#include "headers.h"

extern int **interdiction, **fortification;
extern int Fbudget, Ibudget;
extern int sample, sampleIndex, bindex;
extern double* objSample, M;
extern int checkVisted;
extern int samplingflag;
////////////////////////////////////////////////////////////////////////////////
int Benders(CONFIGURATION param, INSTANCE data, RESULTS *results) {

	//// This section is the Benders code which was available in Github. For the samples, we solve the heuristic algorithm and for the real damage, we use Benders to find the optimal solution //////

	int status = 0;
	int N = data.N;
	int E = data.E;
	int K = data.K;

	CPXENVptr env = NULL;     // pointer to Cplex Enviroment
	CPXLPptr  lp  = NULL;     // pointer to Cplex Linear Program data structure
    
	/////////////////////////////////////////////////////////////////////////////
	// Structs of Global Variables:
	BENDERS_VARIABLES var;
	GLOBAL_INFO	global;
	var.master			= create_double_vector(E+K);
	var.subproblem		= create_double_matrix(K, E+N);
	var.cut_sets		= create_int_matrix(N, N);
	var.capacity		= create_double_matrix(N, N);
	var.dist			= create_double_matrix(N, N);	
	var.last_id				= NONE;
	var.last_id_visits		= 0;
	var.flag_initial_sol	= NO;
	var.flag_generate_cuts	= YES;
	var.flag_fist_master	= YES;
	var.core_point		= create_double_vector(E);
	var.lambda			= 0.0;
	global.data			= &data;
	global.param		= &param;
	global.results		= results;		
	global.var			= &var;
	

	var.init_time = clock();//
	//gettimeofday(&init_time_Linux, NULL);


	if (samplingflag == 0) {
		update_core_point(global, NULL);
	}
	/////////////////////////////////////////////////////////////////////////////
	if (samplingflag == 1) {
		//// For the samples, we solve the heuristic algorithm////
		if (elapsed(var.init_time) < param.MAX_CPU_TIME) {/// /// deactivate this line for calculating time in server
		//if (elapsed(init_time_Linux) < param.MAX_CPU_TIME) { /// activate this line for calculating time in server
			if (param.SCREEN_OUTPUT >= 1) { printf("\nGENERATING BASE MODEL:\n"); }
			status = create_CPLEX_master_lp(env, &lp, global);
			assert(status == 0);
		}
	}

	else {
		//// To solve the OCST problem using benders: 
		// Create the CPLEX Enviroment
		if (elapsed(var.init_time) < param.MAX_CPU_TIME) {
		//if (elapsed(init_time_Linux) < param.MAX_CPU_TIME) {/// activate this line for calculating time in server
			status = create_CPLEX_master_enviroment(&env, param, global);
			assert(status == 0);
		}

		// Generate the initial LP:
		if (elapsed(var.init_time) < param.MAX_CPU_TIME) {
		//if (elapsed(init_time_Linux) < param.MAX_CPU_TIME) {
			if (param.SCREEN_OUTPUT >= 1) { printf("\nGENERATING BASE MODEL:\n"); }
			status = create_CPLEX_master_lp(env, &lp, global);
			assert(status == 0);
		}

		// Now solve the root node:
		if (elapsed(var.init_time) < param.MAX_CPU_TIME) {
		//if (elapsed(init_time_Linux) < param.MAX_CPU_TIME) {
			if (param.SCREEN_OUTPUT >= 1) { printf("\nSOLVING ROOT NODE:\n"); }
			status = solve_root_node(env, lp, global);
			assert(status == 0);

			if (param.SCREEN_OUTPUT >= 1) { printf("\nCLEANING ROOT NODE:\n"); }
			status = clean_root_node(env, lp, global);
			assert(status == 0);

			results->root_time = elapsed(var.init_time);
			//results->root_time = elapsed(init_time_Linux);
		}

		// Solve to integrality
		if (elapsed(var.init_time) < param.MAX_CPU_TIME							// Unless we run out of time...
			&& param.TO_INTEGRALITY == YES										// or we are told to just solve the LP...
			&& results->UpperBound - results->LowerBound > param.MIN_ABS_GAP) {	// or we have already solved the problem...
		//if (elapsed(init_time_Linux) < param.MAX_CPU_TIME							///activate for time caLculation in server. deactivate lines 95-97 
		//	&& param.TO_INTEGRALITY == YES										
		//	&& results->UpperBound - results->LowerBound > param.MIN_ABS_GAP) {
			// Transform the LP into a MILP and set the callbacks
			status = create_CPLEX_master_milp(env, lp, global);
			assert(status == 0);

			if (param.SCREEN_OUTPUT >= 1) { printf("\nSOLVING INTEGER PROBLEM:\n"); }

			// And then solve the problem to integrality
			status = solve_to_integrality(env, lp, global);
			assert(status == 0);
		}

	}
	results->total_time = elapsed(var.init_time);
	//results->total_time = elapsed(init_time_Linux);

	// Clean:
	free_double_vector(&(var.master)); 
	free_double_matrix(&(var.subproblem), K);
	free_int_matrix(&(var.cut_sets), N);
	free_double_matrix(&(var.capacity), N);
	free_double_matrix(&(var.dist), N);
	free_double_vector(&(var.core_point));

	return status;
}
////////////////////////////////////////////////////////////////////////////////
