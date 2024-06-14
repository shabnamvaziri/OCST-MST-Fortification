#include "headers.h"

extern int** interdiction, ** fortification;
extern int Fbudget, Ibudget;
extern int sample, sampleIndex, bindex;
extern int checkVisted;
extern double* ySample, M;
////////////////////////////////////////////////////////////////////////////////
int create_CPLEX_master_enviroment(CPXENVptr* env, CONFIGURATION param, GLOBAL_INFO global) {
	int status = 0;

	// Create enviroment:
	*env = CPXopenCPLEX(&status);
	assert(*env != NULL);

	// Set parameters:
	status = MAX(status, CPXsetintparam(*env, CPX_PARAM_THREADS, 1));						// Threads usados
	status = MAX(status, CPXsetdblparam(*env, CPX_PARAM_TILIM, param.MAX_CPU_TIME));		// Time Limit
	status = MAX(status, CPXsetdblparam(*env, CPX_PARAM_EPGAP, EPSILON * EPSILON));			// Gap de Epsilon Optimalidad
	status = MAX(status, CPXsetdblparam(*env, CPX_PARAM_CUTSFACTOR, 1.0));						// <= 1.0 No cuts will be generated, >1: Limits the number of cuts that can be added
	status = MAX(status, CPXsetintparam(*env, CPX_PARAM_MIPSEARCH, CPX_MIPSEARCH_TRADITIONAL));// Turn on traditional search for use with control callbacks 
	status = MAX(status, CPXsetintparam(*env, CPX_PARAM_MIPCBREDLP, CPX_OFF));					// Let MIP callbacks work on the original model 
	status = MAX(status, CPXsetintparam(*env, CPX_PARAM_HEURFREQ, -1));
	status = MAX(status, CPXsetintparam(*env, CPX_PARAM_PRELINEAR, CPX_OFF));
	status = MAX(status, CPXsetintparam(*env, CPX_PARAM_MIPINTERVAL, 1));
	status = MAX(status, CPXsetintparam(*env, CPX_PARAM_PRESLVND, -1));
	status = MAX(status, CPXsetintparam(*env, CPX_PARAM_PREIND, 0));
	status = MAX(status, CPXsetintparam(*env, CPX_PARAM_REPEATPRESOLVE, 0));
	//status = MAX(status, CPXsetdblparam(*env, CPX_PARAM_CUTUP, (global.results->UpperBound)+0.01));
	status = MAX(status, CPXsetintparam(*env, CPX_PARAM_MIPEMPHASIS, 2));		//	0	CPX_MIPEMPHASIS_BALANCED		Balance optimality and feasibility; default
	//	1	CPX_MIPEMPHASIS_FEASIBILITY		Emphasize feasibility over optimality
	//	2	CPX_MIPEMPHASIS_OPTIMALITY		Emphasize optimality over feasibility
	//	3	CPX_MIPEMPHASIS_BESTBOUND		Emphasize moving best bound
	//	4	CPX_MIPEMPHASIS_HIDDENFEAS		Emphasize finding hidden feasible solutions	
// The Lazy Constraints Callback will switch off these anyway:
	status = MAX(status, CPXsetintparam(*env, CPX_PARAM_REDUCE, CPX_OFF));	// 0: No primal or dual reductions,   1: Only primal reductions ???
	status = MAX(status, CPXsetintparam(*env, CPX_PARAM_PRELINEAR, CPX_OFF));	// Assure linear mappings between the presolved and original models 


	status = MAX(status, CPXsetintparam(*env, CPX_PARAM_VARSEL, 3)); // 0: default, 1: max_infeas, 2: pseudo_cost, 3: strong_branching


	return status;
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
int create_CPLEX_master_lp(CPXENVptr env, CPXLPptr* lp_ptr, GLOBAL_INFO global) {

	clock_t	start;
	int i, j, k, l, c, e, t;
	int added_cuts, status = 0;

	// Just for lazyness:
	int		N = global.data->N;
	int		E = global.data->E;
	int		K = global.data->K;
	int** index_e = global.data->index_e;
	int* index_i = global.data->index_i;
	int* index_j = global.data->index_j;
	int** index_k = global.data->index_k;
	double** W = global.data->W;
	double** C = global.data->C;
	double** C1 = global.data->C1;
	double** C2 = global.data->C2;
	COMMODITY* Com = global.data->Com;
	double longest_path_ub = LongestHamiltonianPathUB(N, C);
	int sssss;
	// New data:
	double** Y = create_double_matrix(N, N);

	// 4-tree data structures:
	double best_tree_value, tree_value;
	double COST[6];
	int TREE_4[16][6] = { {1,1,1,0, 0,0}, {0,1,1,1, 0,0}, {1,0,1,1, 0,0}, {1,1,0,1, 0,0},
						 {1,0,1,0, 1,0}, {0,1,0,1, 0,1}, {1,0,1,0, 0,1}, {0,1,0,1, 1,0},
						 {1,0,0,0, 1,1}, {0,1,0,0, 1,1}, {0,0,1,0, 1,1}, {0,0,0,1, 1,1},
						 {1,1,0,0, 0,1}, {0,1,1,0, 1,0}, {0,0,1,1, 0,1}, {1,0,0,1, 1,0} };

	/// CPLEX MASTER POBLEM ////////////////////////////////////////////////////
	int       numcols; // número de variables ..................................
	int       numrows; // número de restricciones sin contar las cotas .........
	int       numnz;   // número de elementos no nulos de la matriz ............
	int       objsen;  // sentido de la optimizacion (min:1, max:-1 ) ..........
	double* obj;    // coeficientes de la función objetivo ..................
	double* rhs;    // términos independientes de las restricciones .........
	char* sense;  // sentido de las restricciones (<=: 'L', =:'E', >=:'G'). 
	int* matbeg; // índice del primer coeficiente no nulo de cada columna.
	int* matcnt; // número de elementos no nulos de cada columna .........
	int* matind; // fila a la que corresponde cada elemento no nulo  .....
	double* matval; // valores de los coef no nulos en las restricciones ....
	double* lb;     // cotas inferiores de las variables ....................
	double* ub;     // cotas superiores de las variables ....................
	////////////////////////////////////////////////////////////////////////////


	// OPT CUTS FROM HEURISTIC SOLUTIONS
	if (global.param->MST_OPT_CUTS == YES) {		// Minimum Spanning Tree (w.r.t. matrix C)
		MST(N, global.data->C, Y, NULL);
		for (e = 0; e < E; e++) {
			global.var->master[e] = Y[index_i[e]][index_j[e]];
			ySample[(index_i[e] * N + index_j[e]) + (sampleIndex * N * N)] = Y[index_i[e]][index_j[e]];
		}

		free_double_matrix(&Y, N);
		return status;
	}
	////////////////////////////////////////////////////////////////////////////////

}

