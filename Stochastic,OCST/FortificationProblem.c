#include "headers.h"

extern int** interdiction, ** fortification, ** auxiliaryInterdiction, ** waitingInterdiction;
extern int Fbudget, Ibudget, CoverConstraint, sglobal, S;
extern int M, sample, sampleIndex, bindex, iteration, iter, * violation, Fortification_Flag, waitinglist;
extern double* objSample, * ySample, zPerceived;

double FortificationProblem(CONFIGURATION param, INSTANCE data, RESULTS* results) {
	GLOBAL_INFO global;		// Global Paramenters of the algorithm
	global.results = results;
	INSTANCE* instances;		// Array of instance data
	int i, j, k, a, e, numvar, ii;
	int status = 0;
	int N = data.N;
	int E = data.E;
	int K = data.K;
	//double** C = data.C;
	//double** W = data.W;
	//double** Y = create_double_matrix(N, N);
	//double** dist = create_double_matrix(N, N);
	//int** index_e = data.index_e;
	//int* index_i = data.index_i;
	//int* index_j = data.index_j;
	//COMMODITY* Com = data.Com;


	int index, index1, coverindex, waitindex;  // auxiliar indices to fill in the constraint matrix
	double best_upper_bound, best_lower_bound;
	int nodecount;     //Variables to call cplex
	CPXLPptr  lpFP;      // data strucutre to store a problem in cplex ...................
	CPXENVptr envFP;     // cplex environment.............................................
	int       numcols; // number of variables ..........................................
	int       numrows; // number of constraints.........................................
	int       numnz;   // number of non-zero elements in the matrix ....................
	int       objsen;  // optimization sense (min:1, max:-1 ) ..........................
	double* obj;    // objective function coefficients ..............................
	double* rhs;    // right and side of constraints ................................
	char* sense;  // constraints sense (<=: 'L', =:'E', >=:'G') ...................
	int* matbeg; // index of first non-zero element in each row...................
	int* matind; // associated column of each non-zelo element ...................
	double* matval; // coefficient values fo the non-zero elements of constraints....
	double* lb;     // lower bounds of variables.....................................
	double* ub;     // upper bounds of variables.....................................
	double* x;      // solution vector (double, even if the problem is integer) .....
	char      probname[16]; // problem name for cplex .......................................
	char* ctype;  // variable type ('C', 'I', 'B') only if integer.................
	double    value;   // objevtive value of solution ..................................
	double    num_fortification_var;
	int** pos_fortification;

	pos_fortification = create_int_matrix(N, N);

	//Initialize CPLEX environment
	envFP = CPXopenCPLEX(&status);
	if (envFP == NULL) {
		char  errmsg[1024];
		printf("Could not open CPLEX. \n");
		CPXgeterrorstring(envFP, status, errmsg);
		printf("%s", errmsg);
	}

	// Create the problem in CPLEX 
	strcpy(probname, "UFLP");
	lpFP = CPXcreateprob(envFP, &status, probname);
	if (envFP == NULL) {
		char  errmsg[1024];
		printf("Could not create LP. \n");
		CPXgeterrorstring(envFP, status, errmsg);
		printf("%s", errmsg);
	}
	CPXchgobjsen(envFP, lpFP, CPX_MIN);

	///////////////// fortification variable/////////////
	index1 = 0;  // index of columns
	numcols = N * N;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");


	for (i = 0; i < N - 1; i++) {
		for (j = i + 1; j < N; j++) {
			pos_fortification[i][j] = index1;
			obj[index1] = 0;
			ctype[index1] = 'B';
			lb[index1] = 0;
			ub[index1] = 1;
			index1++;
		}
	}
	status = CPXnewcols(envFP, lpFP, index1, obj, lb, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_fortification_var = index1;

	//Add fortification budget constraint 
	numrows = 1;
	numnz = (N * N);
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;

	sense[index1] = 'L';
	rhs[index1] = Fbudget;
	matbeg[index1++] = index;
	for (i = 0; i < N - 1; i++) {
		for (j = i + 1; j < N; j++) {
			matind[index] = pos_fortification[i][j];
			matval[index++] = 1;
		}
	}

	status = CPXaddrows(envFP, lpFP, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	///////////covering constraint/////////////
	numrows = CoverConstraint * S;
	numnz = N * N * CoverConstraint * S;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (sglobal = 0; sglobal < S; sglobal++) {
		for (coverindex = 0; coverindex < CoverConstraint; coverindex++) {
			sense[index1] = 'G';
			rhs[index1] = 1;
			matbeg[index1++] = index;
			for (i = 0; i < N - 1; i++) {
				for (j = i + 1; j < N; j++) {
					matind[index] = pos_fortification[i][j];
					matval[index++] = auxiliaryInterdiction[i * N + j][(coverindex)*S + sglobal];
				}
			}
		}
	}

	status = CPXaddrows(envFP, lpFP, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	///////////covering constraint waiting list/////////////
	numrows = waitinglist * S;
	numnz = N * N * waitinglist * S;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (sglobal = 0; sglobal < S; sglobal++) {
		for (waitindex = 0; waitindex < waitinglist; waitindex++) {
			sense[index1] = 'G';
			rhs[index1] = 1;
			matbeg[index1++] = index;
			for (i = 0; i < N - 1; i++) {
				for (j = i + 1; j < N; j++) {
					matind[index] = pos_fortification[i][j];
					matval[index++] = waitingInterdiction[(i * N + j) + (sglobal * N * N)][waitindex];
				}
			}
		}
	}

	status = CPXaddrows(envFP, lpFP, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);



	CPXwriteprob(envFP, lpFP, "modelFP.lp", NULL);                          //write the model in .lp format if needed (to debug)

	CPXsetintparam(envFP, CPX_PARAM_SCRIND, CPX_OFF); //output display
	CPXsetintparam(envFP, CPX_PARAM_THREADS, 4);
	//CPXsetintparam(envFP,CPX_PARAM_INTSOLLIM,1);    //stops after finding first integer sol.
	CPXsetintparam(envFP, CPX_PARAM_MIPDISPLAY, 3); //different levels of output display
	//CPXsetintparam(envFP,CPX_PARAM_MIPEMPHASIS,1);//0:balanced; 1:feasibility; 2:optimality,3:bestbound, 4:hiddenfeas
	CPXsetdblparam(envFP, CPX_PARAM_TILIM, 86400); // time limit
	//CPXsetdblparam(envFP,CPX_PARAM_TRELIM, 14000); // B&B memory limit
	status = CPXsetintparam(envFP, CPX_PARAM_MEMORYEMPHASIS, 1);	//conserve memory where possible
	CPXsetintparam(envFP, CPX_PARAM_NODEFILEIND, 0);
	CPXsetdblparam(envFP, CPX_PARAM_EPGAP, 0.000001); // e-optimal solution (%gap)
	//CPXsetdblparam(envFP,CPX_PARAM_EPAGAP, 0.0000000001); // e-optimal solution (absolute value)
	//CPXsetdblparam(envFP,CPX_PARAM_EPINT, 0.0000000001); // integer precision
	//CPXsetintparam(envFP,CPX_PARAM_THREADS, 1); // Number of threads to use
	//CPXsetdblparam(envFP,CPX_PARAM_EPRHS, 0.0000001);
	//CPXsetintparam(envFP,CPX_PARAM_REDUCE, 0);  // only needed when adding lazy constraints
	//CPXsetintparam(envFP,CPX_PARAM_HEURFREQ, -1); //heuristic frequency and intensisty 
	//CPXsetdblparam(envFP, CPX_PARAM_CUTSFACTOR, 1.0);  //limit the number of cuts added by cplex 1.0002
	//CPXsetdblparam(envFP,CPX_PARAM_CUTUP,UpperBound+.01); // provide an initial upper bound
	//CPXsetintparam(envFP,CPX_PARAM_MIPEMPHASIS,CPX_MIPEMPHASIS_OPTIMALITY);  // MIP emphasis: optimality, feasibility, moving best bound
	//CPXsetintparam(envFP,CPX_PARAM_PARALLELMODE, 1); 
	//CPXsetintparam(envFP,CPX_PARAM_PREIND,0);
	//CPXsetintparam(envFP,CPX_PARAM_MIPORDIND,CPX_ON); // Turn on or off the use of priorities on bracnhing variables
	//CPXsetintparam(envFP,CPX_PARAM_MIPEMPHASIS,CPX_MIPEMPHASIS_BESTBOUND);  // MIP emphasis: optimality, feasibility, moving best bound

	//gettimeofday(&start, NULL);
	//startTemp = clock();
	CPXmipopt(envFP, lpFP);  //solve the integer program
	//gettimeofday(&stop, NULL);
	//endTemp = clock();
	//MasterTime = (double)(endTemp - startTemp) / (double)(CLOCKS_PER_SEC);
	//MasterTime = ((double)(stop.tv_sec - start.tv_sec) * 1000 + (double)(stop.tv_usec - start.tv_usec) / 1000) / 1000;

	i = CPXgetstat(envFP, lpFP);
	if (i == 101)
		printf("Optimal solution found\n");
	else if (i == 102)
		printf("e-optimal solution found\n");
	else if (i == 103) {
		printf(" infeasible solution\n");
		Fortification_Flag = 1;
	}
	else if (i == 107)
		printf("Time limit reached\n");
	else
		printf("Unknown stopping criterion (%d)\n", i);

	// retrive solution values
	CPXgetmipobjval(envFP, lpFP, &value);
	printf("Upper bound: %.2f   ", value);
	best_upper_bound = value;
	//best_upper_bound += 2 * x1F + x2F;
	// If CPLEX was able to find the optimal solution, the previous function provides the optimal solution value
	//if not, it provides the best upper bound
	CPXgetbestobjval(envFP, lpFP, &value);  //best lower bound in case the problem was not solved to optimality
	best_lower_bound = value;
	printf("Lower bound: %.2f  \n", value);

	nodecount = CPXgetnodecnt(envFP, lpFP);
	printf("Number of BB nodes : %ld  \n", nodecount);

	numcols = CPXgetnumcols(envFP, lpFP);
	d_vector(&x, numcols, "open_cplex:0");
	CPXgetmipx(envFP, lpFP, x, 0, numcols - 1);  // obtain the values of the decision variables

	if (lpFP != NULL) {
		status = CPXfreeprob(envFP, &lpFP);
		if (status) {
			fprintf(stderr, "CPXfreeprob failed, error code %d.\n", status);
		}
	}
	if (envFP != NULL) {
		status = CPXcloseCPLEX(&envFP);
		if (status) {
			char  errmsg[1024];
			fprintf(stderr, "Could not close CPLEX envFPironment.\n");
			CPXgeterrorstring(envFP, status, errmsg);
			fprintf(stderr, "%s", errmsg);
		}
	}

	index = 0;

	for (i = 0; i < N - 1; i++) {
		for (j = i + 1; j < N; j++) {
			//	printf("i= %d, \t j=%d \t index=%d\n", i, j, index);

			if (x[index] > 0.5) fortification[i][j] = 1;
			if (x[index] <= 0.5) fortification[i][j] = 0;
			index++;
			if (fortification[i][j] > 0) {
				printf("fortification[%d][%d]=%d\n", i, j, fortification[i][j]);
			}
		}
	}

	free(x);
	free(pos_fortification);
	return best_upper_bound;
}
