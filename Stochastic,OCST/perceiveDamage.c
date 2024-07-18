#include "headers.h"

extern int** interdiction, ** fortification, ** auxiliaryInterdiction, ** interdictionperceived, ** InterdictionS;
extern int Fbudget, *Ibudget, CoverConstraint, sglobal;
extern int sample, sampleIndex, bindex, iteration, iter, * violation;
extern double* objSample, * ySample, zPerceived, M, Ubound;

double perceivedDamage(CONFIGURATION param, INSTANCE data, RESULTS* results) {
	GLOBAL_INFO global;		// Global Paramenters of the algorithm
	global.results = results;
	INSTANCE* instances;		// Array of instance data
	int i, j, k, a, e, numvar, ii, indexsample;
	int status = 0;
	int N = data.N;
	int E = data.E;
	int K = data.K;
	double** C = data.C;
	//double** W = data.W;
	//double** Y = create_double_matrix(N, N);
	//double** dist = create_double_matrix(N, N);
	int** index_e = data.index_e;
	int* index_i = data.index_i;
	int* index_j = data.index_j;
	COMMODITY* Com = data.Com;


	int index, index1;  // auxiliar indices to fill in the constraint matrix
	double best_upper_bound, best_lower_bound;
	int nodecount;     //Variables to call cplex
	CPXLPptr  lpPD;      // data strucutre to store a problem in cplex ...................
	CPXENVptr envPD;     // cplex environment.............................................
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
	double    num_z_var, num_interdiction_var;
	double    pos_z;
	int** pos_interdiction;

	pos_interdiction = create_int_matrix(N, N);

	//Initialize CPLEX environment
	envPD = CPXopenCPLEX(&status);
	if (envPD == NULL) {
		char  errmsg[1024];
		printf("Could not open CPLEX. \n");
		CPXgeterrorstring(envPD, status, errmsg);
		printf("%s", errmsg);
	}

	// Create the problem in CPLEX 
	strcpy(probname, "UFLP");
	lpPD = CPXcreateprob(envPD, &status, probname);
	if (envPD == NULL) {
		char  errmsg[1024];
		printf("Could not create LP. \n");
		CPXgeterrorstring(envPD, status, errmsg);
		printf("%s", errmsg);
	}
	CPXchgobjsen(envPD, lpPD, CPX_MAX);

	//Define z variables
	index1 = 0;  // index of columns
	numcols = 1;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");

	pos_z = index1;
	obj[index1] = 1;
	ctype[index1] = 'C';
	lb[index1] = objSample[0];
	ub[index1] = objSample[0] + (Ibudget[sglobal]) * M;
	index1++;

	status = CPXnewcols(envPD, lpPD, index1, obj, lb, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_z_var = index1;

	///////////////// interdiction variable/////////////
	index1 = 0;  // index of columns
	numcols = N * N;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");

	for (i = 0; i < N - 1; i++) {
		for (j = i + 1; j < N; j++) {
			pos_interdiction[i][j] = index1 + num_z_var;
			obj[index1] = 0;
			ctype[index1] = 'B';
			lb[index1] = 0;
			ub[index1] = 1;
			index1++;
		}
	}
	status = CPXnewcols(envPD, lpPD, index1, obj, lb, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_interdiction_var = index1;

	//Add constraint 1 // z <= cost(sample)+ M*x
	numrows = sampleIndex;
	numnz = (1 + N * N) * sampleIndex * 2;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (indexsample = 0; indexsample < sampleIndex; indexsample++) {
		sense[index1] = 'L';
		rhs[index1] = objSample[indexsample];
		matbeg[index1++] = index;
		matind[index] = pos_z;
		matval[index++] = 1;
		for (i = 0; i < N - 1; i++) {
			for (j = i + 1; j < N; j++) {
				if (ySample[(i * N + j) + (indexsample * N * N)] > 0.5) {
					matind[index] = pos_interdiction[i][j];
					matval[index++] = -M;
				}
			}
		}
	}

	status = CPXaddrows(envPD, lpPD, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	///////////Interdiction budget/////////////
	numrows = 1 * N;
	numnz = N * N;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;

	sense[index1] = 'L';
	rhs[index1] = Ibudget[sglobal];
	matbeg[index1++] = index;
	for (i = 0; i < N - 1; i++) {
		for (j = i + 1; j < N; j++) {
			matind[index] = pos_interdiction[i][j];
			matval[index++] = 1;
		}
	}

	status = CPXaddrows(envPD, lpPD, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	///////////interdiction < 1 - fortification
	numrows = N * N;
	numnz = (N * N) * N * N;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (i = 0; i < N - 1; i++) {
		for (j = i + 1; j < N; j++) {
			sense[index1] = 'L';
			rhs[index1] = 1 - fortification[i][j];
			matbeg[index1++] = index;
			matind[index] = pos_interdiction[i][j];
			matval[index++] = 1;
		}
	}

	status = CPXaddrows(envPD, lpPD, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);


	//CPXwriteprob(envPD, lpPD, "model2.lp", NULL);                          //write the model in .lp format if needed (to debug)

	CPXsetintparam(envPD, CPX_PARAM_SCRIND, CPX_OFF); //output display
	CPXsetintparam(envPD, CPX_PARAM_THREADS, 4);
	//CPXsetintparam(envPD,CPX_PARAM_INTSOLLIM,1);    //stops after finding first integer sol.
	CPXsetintparam(envPD, CPX_PARAM_MIPDISPLAY, 3); //different levels of output display
	//CPXsetintparam(envPD,CPX_PARAM_MIPEMPHASIS,1);//0:balanced; 1:feasibility; 2:optimality,3:bestbound, 4:hiddenfeas
	CPXsetdblparam(envPD, CPX_PARAM_TILIM, 86400); // time limit
	//CPXsetdblparam(envPD,CPX_PARAM_TRELIM, 14000); // B&B memory limit
	status = CPXsetintparam(envPD, CPX_PARAM_MEMORYEMPHASIS, 1);	//conserve memory where possible
	CPXsetintparam(envPD, CPX_PARAM_NODEFILEIND, 0);
	CPXsetdblparam(envPD, CPX_PARAM_EPGAP, 0.000001); // e-optimal solution (%gap)
	//CPXsetdblparam(envPD,CPX_PARAM_EPAGAP, 0.0000000001); // e-optimal solution (absolute value)
	//CPXsetdblparam(envPD,CPX_PARAM_EPINT, 0.0000000001); // integer precision
	//CPXsetintparam(envPD,CPX_PARAM_THREADS, 1); // Number of threads to use
	//CPXsetdblparam(envPD,CPX_PARAM_EPRHS, 0.0000001);
	//CPXsetintparam(envPD,CPX_PARAM_REDUCE, 0);  // only needed when adding lazy constraints
	//CPXsetintparam(envPD,CPX_PARAM_HEURFREQ, -1); //heuristic frequency and intensisty 
	//CPXsetdblparam(envPD, CPX_PARAM_CUTSFACTOR, 1.0);  //limit the number of cuts added by cplex 1.0002
	//CPXsetdblparam(envPD,CPX_PARAM_CUTUP,UpperBound+.01); // provide an initial upper bound
	//CPXsetintparam(envPD,CPX_PARAM_MIPEMPHASIS,CPX_MIPEMPHASIS_OPTIMALITY);  // MIP emphasis: optimality, feasibility, moving best bound
	//CPXsetintparam(envPD,CPX_PARAM_PARALLELMODE, 1); 
	//CPXsetintparam(envPD,CPX_PARAM_PREIND,0);
	//CPXsetintparam(envPD,CPX_PARAM_MIPORDIND,CPX_ON); // Turn on or off the use of priorities on bracnhing variables
	//CPXsetintparam(envPD,CPX_PARAM_MIPEMPHASIS,CPX_MIPEMPHASIS_BESTBOUND);  // MIP emphasis: optimality, feasibility, moving best bound

	//gettimeofday(&start, NULL);
	//startTemp = clock();
	CPXmipopt(envPD, lpPD);  //solve the integer program
	//gettimeofday(&stop, NULL);
	//endTemp = clock();
	//MasterTime = (double)(endTemp - startTemp) / (double)(CLOCKS_PER_SEC);
	//MasterTime = ((double)(stop.tv_sec - start.tv_sec) * 1000 + (double)(stop.tv_usec - start.tv_usec) / 1000) / 1000;

	i = CPXgetstat(envPD, lpPD);
	if (i == 101)
		printf("Optimal solution found\n");
	else if (i == 102)
		printf("e-optimal solution found\n");
	else if (i == 103)
		printf(" infeasible solution\n");
	else if (i == 107)
		printf("Time limit reached\n");
	else
		printf("Unknown stopping criterion (%d)\n", i);

	// retrive solution values
	CPXgetmipobjval(envPD, lpPD, &value);
	printf("Upper bound: %.2f   ", value);
	best_upper_bound = value;
	//best_upper_bound += 2 * x1F + x2F;
	// If CPLEX was able to find the optimal solution, the previous function provides the optimal solution value
	//if not, it provides the best upper bound
	CPXgetbestobjval(envPD, lpPD, &value);  //best lower bound in case the problem was not solved to optimality
	best_lower_bound = value;
	printf("Lower bound: %.2f  \n", value);

	nodecount = CPXgetnodecnt(envPD, lpPD);
	printf("Number of BB nodes : %ld  \n", nodecount);

	numcols = CPXgetnumcols(envPD, lpPD);
	d_vector(&x, numcols, "open_cplex:0");
	CPXgetmipx(envPD, lpPD, x, 0, numcols - 1);  // obtain the values of the decision variables

	if (lpPD != NULL) {
		status = CPXfreeprob(envPD, &lpPD);
		if (status) {
			fprintf(stderr, "CPXfreeprob failed, error code %d.\n", status);
		}
	}
	if (envPD != NULL) {
		status = CPXcloseCPLEX(&envPD);
		if (status) {
			char  errmsg[1024];
			fprintf(stderr, "Could not close CPLEX envPDironment.\n");
			CPXgeterrorstring(envPD, status, errmsg);
			fprintf(stderr, "%s", errmsg);
		}
	}

	index = 0;

	Ubound = x[index];
	index++;

	for (i = 0; i < N - 1; i++) {
		for (j = i + 1; j < N; j++) {
			if (x[index] > 0.5) interdiction[i][j] = 1;
			if (x[index] <= 0.5) interdiction[i][j] = 0;
			interdictionperceived[i][j] = interdiction[i][j];
			index++;
			if (interdiction[i][j] > 0) {
				printf("index=%d\t", index - 1);
				printf("interdiction[%d][%d]=%d\n", i, j, interdiction[i][j]);
			}
		}
	}

	for (i = 0; i < N - 1; i++) {
		for (j = i + 1; j < N; j++) {
			InterdictionS[i * N + j][sglobal] = interdiction[i][j];
		}
	}
	free(x);
	free(pos_interdiction);

	return best_upper_bound;
}
