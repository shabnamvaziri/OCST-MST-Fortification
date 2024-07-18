#include "headers.h"

extern int** interdiction, ** fortification, ** auxiliaryInterdiction, ** incumbantInterdiction, ** incumbantfortification, ** incumbantY, ** yAuxiliary;
extern int Fbudget, *Ibudget, Fortification_Flag, CoverConstraint, S, sglobal;
extern int M, sample, sampleIndex, bindex, iteration, iter, sampleRemoved, acceptedSample, FlagCover;
extern double* objSample, * ySample, zPerceived, lowerboundsto, upperboundsto, zbar;
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
	//double** C = data.C;
	//double** W = data.W;
	//double** Y = create_double_matrix(N, N);
	//double** dist = create_double_matrix(N, N);
	//int** index_e = data.index_e;
	//int* index_i = data.index_i;
	//int* index_j = data.index_j;
	//COMMODITY* Com = data.Com;


	Out = open_file(outfile, "a+");
	fprintf(Out, "%s\t", data.filename);
	fprintf(Out, "N: %d\t", data.N);
	fprintf(Out, "E: %d\t", data.E);
	fprintf(Out, "R: %d\t", data.K);
	fprintf(Out, "Fbudget: %d\t", Fbudget);
	fprintf(Out, "S: %d\t", S);
	fprintf(Out, "sample: %d\t", sample);

	for (sglobal = 0; sglobal < S; sglobal++) {
		for (i = 0; i < N - 1; i++) {
			for (j = i + 1; j < N; j++) {
				if (incumbantInterdiction[i * N + j][sglobal] > 0) {
					fprintf(Out, "interdiction[%d][%d][%d]\t", i, j, sglobal);
				}
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
	for (sglobal = 0; sglobal < S; sglobal++) {
		for (i = 0; i < N - 1; i++) {
			for (j = i + 1; j < N; j++) {
				if (incumbantY[i * N + j][sglobal] > 0) {
					fprintf(Out, "y[%d][%d][%d]\t", i, j, sglobal);
				}
			}
		}
	}

	fprintf(Out, "zbar = %f \t", incumbantzbar);
	fprintf(Out, "UB = %f \t", incumbantUB);
	fprintf(Out, "LB = %f \t", incumbantLB);
	fprintf(Out, "zfinal = %f \t", zbar);
	fprintf(Out, "UBfinal = %f \t", upperboundsto);
	fprintf(Out, "LBfinal = %f \t", lowerboundsto);
	fprintf(Out, "iteration = %d \t", sampleIndex);
	fprintf(Out, "fortificationFlag = %d \t", Fortification_Flag);


	fclose(Out);


	return 0;
}

