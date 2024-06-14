#include "headers.h"

extern int** interdiction, ** fortification, ** yAuxiliary;
extern int Fbudget, Ibudget;
extern int sample, sampleIndex, bindex, iteration, iter;
extern double* objSample, * ySample, zPerceived, UBound, LBound, M;

double realDamage(CONFIGURATION param, INSTANCE data, RESULTS* results) {
	GLOBAL_INFO global;		// Global Paramenters of the algorithm
	BENDERS_VARIABLES var;
	global.var = &var;
	global.results = results;
	INSTANCE* instances;		// Array of instance data
	int i, j, k, a, e, numvar, ii;
	int status = 0;
	int N = data.N;
	int E = data.E;
	int K = data.K;
	double** C = data.C;
	double** W = data.W;
	
	int** index_e = data.index_e;
	int* index_i = data.index_i;
	int* index_j = data.index_j;
	COMMODITY* Com = data.Com;

	if (param.ALGORITHM == 0 || param.ALGORITHM == 1) {		// Using Benders
		status = Benders(param, data, results);
		assert(status == 0);
	}
	for (i = 0; i < N - 1; i++) {             // Y_ij variables...
		for (j = i + 1; j < N; j++) {
			objSample[sampleIndex] += ySample[(i * N + j) + (sampleIndex * N * N)] * data.C[i][j];
		}
	}

	LBound = objSample[sampleIndex];
	for (e = 0; e < E; e++) {
		yAuxiliary[index_i[e]][index_j[e]] = ySample[(index_i[e] * N + index_j[e]) + (sampleIndex * N * N)];
	}
}


