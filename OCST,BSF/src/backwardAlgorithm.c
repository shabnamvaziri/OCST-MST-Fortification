#include "headers.h"

extern int** interdiction, ** fortification, ** auxiliaryInterdiction, **incumbantInterdiction, ** incumbantfortification, **incumbantY, **yAuxiliary;
extern int Fbudget, Ibudget, Fortification_Flag, CoverConstraint;
extern int sample, sampleIndex, bindex, iteration, iter, sampleRemoved, acceptedSample, FlagCover;
extern double* objSample, * ySample, zPerceived, UBound, LBound, zbar;
extern int checkVisted;
extern int iglobal, initialsample;
extern double incumbantzbar, incumbantUB, incumbantLB, M;
clock_t			  startTemp, endTemp;

extern FILE* Out;
extern char OutName[100];
extern char outfile[20];


double backwardAlgorithm(CONFIGURATION param, INSTANCE data, RESULTS* results) {
	GLOBAL_INFO global;		// Global Paramenters of the algorithm
	global.results = results;
	INSTANCE* instances;		// Array of instance data
	int i, j, k, a, e, numvar, ii, index1;
	int status = 0;
	int N = data.N;
	int E = data.E;
	int K = data.K;
	


	for (i = 0; i < N - 1; i++) {
		for (j = i + 1; j < N; j++) {
			incumbantInterdiction[i][j] = auxiliaryInterdiction[i * N + j][CoverConstraint - 1];
			//printf("interdiction[%d][%d]=%d\n", i, j, incumbantInterdiction[i][j]);
		}
	}

	for (i = 0; i < N - 1; i++) {
		for (j = i + 1; j < N; j++) {
			incumbantfortification[i][j] = fortification[i][j];
		}
	}

	for (i = 0; i < N - 1; i++) {
		for (j = i + 1; j < N; j++) {
			incumbantY[i][j] = yAuxiliary[i][j];
		}
	}
	incumbantzbar = zbar;
	incumbantUB = UBound;
	incumbantLB = LBound;

	return 0;
}

