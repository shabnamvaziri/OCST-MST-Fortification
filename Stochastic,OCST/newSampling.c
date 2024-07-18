#include "headers.h"

extern int** interdiction, ** fortification;
extern int Fbudget, *Ibudget;
extern int sample, sampleIndex, bindex, iteration, iter, sample_index;
extern double* objSample, * ySample, M;

int newSampling(CONFIGURATION param, INSTANCE data, RESULTS* results) {
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
	/*double** W = data.W;
	double** Y = create_double_matrix(N, N);
	double** dist = create_double_matrix(N, N);*/
	int** index_e = data.index_e;
	int* index_i = data.index_i;
	int* index_j = data.index_j;
	COMMODITY* Com = data.Com;

	int bsample, bRandom, * bSelected, bFlag;

	bSelected = create_int_vector(sample);
	bsample = 0;

	bFlag = 0;
	for (e = 0; e < E; e++) {
		interdiction[index_i[e]][index_j[e]] = 0;
	}
	bsample = (N - 2);
	for (bindex = 0; bindex < bsample; bindex++) {
		bFlag = 0;
		bRandom = rand() % E;
		bSelected[bindex] = bRandom;
		for (ii = 0; ii < bindex; ii++) {
			if (bSelected[ii] == bRandom) {
				bFlag = 1;
				break;
			}
		}
		if (bFlag == 1) {
			bindex--;
		}
		else {
			interdiction[index_i[bRandom]][index_j[bRandom]] = 1;
		}
	}

	if (param.ALGORITHM == 0 || param.ALGORITHM == 1) {		// Using Benders
		status = Benders(param, data, results);
		assert(status == 0);
	}
	for (i = 0; i < N - 1; i++) {             // Y_ij variables...
		for (j = i + 1; j < N; j++) {
			objSample[sample_index] += ySample[(i * N + j) + (sample_index * N * N)] * data.C[i][j];
		}
	}
	objSample[sample_index] = global.results->UpperBound;
	for (e = 0; e < E; e++) {
		ySample[(index_i[e] * N + index_j[e]) + (sample_index * N * N)] = global.results->best_solution[e];
	}

	return bFlag;
}