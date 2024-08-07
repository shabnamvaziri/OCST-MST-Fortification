#include "headers.h"

extern int** interdiction, ** fortification, ** auxiliaryInterdiction, ** interdictionperceived, ** waitingInterdiction;
extern int Fbudget, Ibudget, CoverConstraint, FlagCover, indexGlobal;
extern int M, sample, sampleIndex, bindex, iteration, iter, * violation, Fortification_Flag, sglobal, S;
extern double* objSample, * ySample, zPerceived, UBound, LBound;

double addCoverconstraint(CONFIGURATION param, INSTANCE data, RESULTS* results) {
	GLOBAL_INFO global;		// Global Paramenters of the algorithm
	global.results = results;
	INSTANCE* instances;		// Array of instance data
	int i, j, k, a, e, numvar, ii;
	int status = 0;
	int N = data.N;
	int E = data.E;
	int K = data.K;
	/*double** C = data.C;
	double** W = data.W;
	double** Y = create_double_matrix(N, N);
	double** dist = create_double_matrix(N, N);
	int** index_e = data.index_e;
	int* index_i = data.index_i;
	int* index_j = data.index_j;
	COMMODITY* Com = data.Com;*/

	for (sglobal = 0; sglobal < S; sglobal++) {
		for (i = 0; i < N - 1; i++) {
			for (j = i + 1; j < N; j++) {
				auxiliaryInterdiction[(i * N + j)][CoverConstraint * S + sglobal] = waitingInterdiction[(i * N + j) + (sglobal * N * N)][indexGlobal];
			}
		}
	}
	CoverConstraint++;
	FlagCover = 1;

	return 0;
}
