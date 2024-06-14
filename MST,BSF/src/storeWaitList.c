#include "headers.h"

extern int** interdiction, ** fortification, ** auxiliaryInterdiction, ** incumbantInterdiction, ** incumbantfortification, ** incumbantY, ** yAuxiliary, ** interdictionperceived;
extern int Fbudget, Ibudget, Fortification_Flag, CoverConstraint, ** waitingInterdiction, waitinglist, ** waitingfortification;
extern int M, sample, sampleIndex, bindex, iteration, iter, sampleRemoved, acceptedSample, FlagCover;
extern double* objSample, * ySample, zPerceived, UBound, LBound, zbar, * waitinglistRealdamage;
extern int checkVisted;
extern int iglobal, initialsample;
extern double incumbantzbar, incumbantUB, incumbantLB;
clock_t			  startTemp, endTemp;

extern FILE* Out;
extern char OutName[100];
extern char outfile[20];


double storeWaitList(CONFIGURATION param, INSTANCE data, RESULTS* results) {
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
			waitingInterdiction[i * N + j][waitinglist] = interdictionperceived[i][j];
			waitingfortification[i * N + j][waitinglist] = fortification[i][j];
		}
	}
	waitinglistRealdamage[waitinglist] = LBound;


	///  I still need to understand why to store w


	waitinglist++;

	return 0;
}

