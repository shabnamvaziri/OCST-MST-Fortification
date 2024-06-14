#include "headers.h"

extern int** interdiction, ** fortification, ** auxiliaryInterdiction, ** incumbantInterdiction, ** incumbantfortification, ** incumbantY, ** yAuxiliary;
extern int Fbudget, Ibudget, Fortification_Flag, CoverConstraint;
extern int M, sample, sampleIndex, bindex, iteration, iter, sampleRemoved, acceptedSample, FlagCover;
extern double* objSample, * ySample, zPerceived, UBound, LBound, zbar;
extern int checkVisted;
extern int iglobal, initialsample;
extern double incumbantzbar, incumbantUB, incumbantLB;
clock_t			  startTemp, endTemp;

extern FILE* Out;
extern char OutName[100];
extern char outfile[20];


double backwardAlgorithmfinal(CONFIGURATION param, INSTANCE data, RESULTS* results) {
	GLOBAL_INFO global;		// Global Paramenters of the algorithm
	global.results = results;
	INSTANCE* instances;		// Array of instance data
	int i, j, k, a, e, numvar, ii, index1;
	int status = 0;
	int N = data.N;
	int E = data.E;
	int K = data.K;
	
	Out = open_file(outfile, "a+");
	fprintf(Out, "%s\t", data.filename);
	fprintf(Out, "N: %d\t", data.N);
	fprintf(Out, "E: %d\t", data.E);
	fprintf(Out, "R: %d\t", data.K);
	fprintf(Out, "Fbudget: %d\t", Fbudget);
	fprintf(Out, "Ibudget: %d\t", Ibudget);
	fprintf(Out, "sample: %d\t", initialsample);
	for (i = 0; i < N - 1; i++) {
		for (j = i + 1; j < N; j++) {
			if (incumbantInterdiction[i][j] > 0) {
				fprintf(Out, "interdiction[%d][%d]\t", i, j);
			}
		}
	}

	for (i = 0; i < N - 1; i++) {
		for (j = i + 1; j < N; j++) {
			if (incumbantfortification[i][j] > 0) {
				fprintf(Out, "fortification[%d][%d]\t", i, j);
			}
		}
	}

	for (i = 0; i < N - 1; i++) {
		for (j = i + 1; j < N; j++) {
			if (incumbantY[i][j] > 0) {
				fprintf(Out, "y[%d][%d]\t", i, j);
			}
		}
	}

	fprintf(Out, "zbar = %f \t", incumbantzbar);
	fprintf(Out, "UB = %f \t", incumbantUB);
	fprintf(Out, "LB = %f \t", incumbantLB);
	fprintf(Out, "zfinal = %f \t", zbar);
	fprintf(Out, "UBfinal = %f \t", UBound);
	fprintf(Out, "LBfinal = %f \t", LBound);
	fprintf(Out, "iteration = %d \t", sampleIndex);
	fprintf(Out, "fortificationFlag = %d \t", Fortification_Flag);


	fclose(Out);


	return 0;
}

