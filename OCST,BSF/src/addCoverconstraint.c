#include "headers.h"

extern int** interdiction, ** fortification, ** auxiliaryInterdiction, ** interdictionperceived, ** waitingInterdiction;
extern int Fbudget, Ibudget, CoverConstraint, FlagCover, indexGlobal;
extern int M, sample, sampleIndex, bindex, iteration, iter, * violation, Fortification_Flag;
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


	for (i = 0; i < N - 1; i++) {
		for (j = i + 1; j < N; j++) {
			auxiliaryInterdiction[i * N + j][CoverConstraint] = waitingInterdiction[i * N + j][indexGlobal];
		}
	}
	CoverConstraint++;
	FlagCover = 1;

	return 0;
}
