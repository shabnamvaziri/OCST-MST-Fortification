#include "headers.h"

extern int** interdiction, ** fortification, ** auxiliaryInterdiction, ** incumbantInterdiction, ** incumbantfortification, ** incumbantY, ** yAuxiliary;
extern int Fbudget, *Ibudget, Fortification_Flag, CoverConstraint, sglobal, S;
extern int sample, sampleIndex, bindex, iteration, iter, sampleRemoved, acceptedSample, FlagCover;
extern double* objSample, * ySample, zPerceived, lowerboundsto, upperboundsto, zbar;
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
				incumbantInterdiction[i * N + j][sglobal] = auxiliaryInterdiction[i * N + j][(CoverConstraint - 1) * S + sglobal];
				//printf("interdiction[%d][%d]=%d\n", i, j, incumbantInterdiction[i][j]);
			}
		}
	}

	for (i = 0; i < N - 1; i++) {
		for (j = i + 1; j < N; j++) {
			incumbantfortification[i][j] = fortification[i][j];
		}
	}
	for (sglobal = 0; sglobal < S; sglobal++) {
		for (i = 0; i < N - 1; i++) {
			for (j = i + 1; j < N; j++) {
				incumbantY[i * N + j][sglobal] = yAuxiliary[i * N + j][sglobal];
			}
		}
	}
	incumbantzbar = zbar;
	incumbantUB = upperboundsto;
	incumbantLB = lowerboundsto;
	//free_double_matrix(&dist, N);
	//free_double_matrix(&Y, N);
	return 0;
}

