#include "headers.h"

extern int** interdiction, ** fortification, ** auxiliaryInterdiction, ** incumbantInterdiction, ** incumbantfortification, ** incumbantY, ** yAuxiliary, ** interdictionperceived;
extern int Fbudget, Ibudget, Fortification_Flag, CoverConstraint, ** waitingInterdiction, waitinglist, ** waitingfortification, ** InterdictionS;
extern int M, sample, sampleIndex, bindex, iteration, iter, sampleRemoved, acceptedSample, FlagCover, sglobal, S;
extern double* objSample, * ySample, zPerceived, UBound, lowerboundsto, zbar, * waitinglistRealdamage;
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
	//double** C = data.C;
	//double** W = data.W;
	//double** Y = create_double_matrix(N, N);
	//double** dist = create_double_matrix(N, N);
	//int** index_e = data.index_e;
	//int* index_i = data.index_i;
	//int* index_j = data.index_j;
	//COMMODITY* Com = data.Com;

	for (sglobal = 0; sglobal < S; sglobal++) {
		for (i = 0; i < N - 1; i++) {
			for (j = i + 1; j < N; j++) {
				waitingInterdiction[(i * N + j) + (sglobal * N * N)][waitinglist] = InterdictionS[i * N + j][sglobal];
				waitingfortification[(i * N + j) + (sglobal * N * N)][waitinglist] = fortification[i][j];
			}
		}
	}
	waitinglistRealdamage[waitinglist] = lowerboundsto;

	waitinglist++;

	return 0;
}

