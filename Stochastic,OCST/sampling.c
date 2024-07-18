#include "headers.h"

extern int** interdiction, ** fortification;
extern int Fbudget, *Ibudget;
extern int sample, sampleIndex, bindex, iteration, iter;
extern double* objSample, * ySample, M;
extern int checkVisted;
int* flag;
extern int** countY, MaxCountY;

extern double TotalTime, SampleTimeLimit, checktime;
extern FILE* Out;
extern char OutName[100];
extern char outfile[20];


int sampling(CONFIGURATION param, INSTANCE data, RESULTS* results) {
	GLOBAL_INFO global;			 // Global Paramenters of the algorithm
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
	checktime = 0;


	MaxCountY = (int)30;
	countY = create_int_matrix(N, N);
	for (i = 0; i < N - 1; i++) {
		for (j = i + 1; j < N; j++) {
			countY[i][j] = 0;
			countY[j][i] = 0;
		}
	}

	COMMODITY* Com = data.Com;
	sampleIndex = 0;
	int bsample, bRandom, * bSelected, bFlag;
	//double stopsample, startsample;
	//checktime = 0;

	bSelected = create_int_vector(sample);
	bsample = 0;
	
	if (param.ALGORITHM == 0 || param.ALGORITHM == 1) {		// Using Benders
		status = Benders(param, data, results);
		assert(status == 0);
	}

	for (i = 0; i < N - 1; i++) {             // Y_ij variables...
		for (j = i + 1; j < N; j++) {
			if (ySample[(i * N + j) + (sampleIndex * N * N)] > 0) {
				objSample[sampleIndex] += ySample[(i * N + j) + (sampleIndex * N * N)] * data.C[i][j];
				countY[i][j]++;
				countY[j][i]++;
			}
		}
	}
	Sorted = (CostSort*)calloc(N * N, sizeof(CostSort));

	int counter = 0;
	for (i = 0; i < N - 1; i++) {             // Y_ij variables...
		for (j = i + 1; j < N; j++) {
			if (ySample[(i * N + j) + (0 * N * N)] == 1) {
				Sorted[counter].index = index_e[i][j];
				Sorted[counter].value = data.C[i][j];
				counter++;
			}
		}
	}
	qsort((CostSort*)Sorted, counter, sizeof(Sorted[0]), Comparevalue);

	sampleIndex++;

	gettimeofday(&startsample, NULL);

	//startsample = clock();//

	for (sampleIndex = 1; sampleIndex < sample; sampleIndex++) {

		bFlag = 0;
		for (e = 0; e < E; e++) {
			interdiction[index_i[e]][index_j[e]] = 0;
		}

		bsample = (rand() % (N - 2) + 1);
		//bsample = fmin(rand() % (N - 2) + 1,Ibudget);

		for (bindex = 0; bindex < bsample; bindex++) {
			interdiction[index_i[Sorted[bindex].index]][index_j[Sorted[bindex].index]] = 1;
			printf("interdiction[%d][%d] = %d\n", index_i[Sorted[bindex].index], index_j[Sorted[bindex].index], interdiction[index_i[Sorted[bindex].index]][index_j[Sorted[bindex].index]]);
		}

		/*status = Heuristics(param, data, results);
		assert(status == 0);*/

		if (param.ALGORITHM == 0 || param.ALGORITHM == 1) {		// Using Benders
			status = Benders(param, data, results);
			assert(status == 0);
		}

		gettimeofday(&stopsample, NULL);
		checktime = ((double)(stopsample.tv_sec - startsample.tv_sec) * 1000 + (double)(stopsample.tv_usec - startsample.tv_usec) / 1000) / 1000;

		//stopsample = clock();//
		//checktime = (double)(stopsample - startsample) / (double)(CLOCKS_PER_SEC);//

		for (i = 0; i < N - 1; i++) {             // Y_ij variables...
			for (j = i + 1; j < N; j++) {
				objSample[sampleIndex] += ySample[(i * N + j) + (sampleIndex * N * N)] * data.C[i][j];
			}
		}
		counter = 0;
		for (e = 0; e < E; e++) {
			if (ySample[(index_i[e] * N + index_j[e]) + (sampleIndex * N * N)] == 1) {
				Sorted[counter].index = e;
				Sorted[counter].value = data.C[index_i[e]][index_j[e]];
				counter++;
			}
		}
		for (e = 0; e < E; e++) {
			if (ySample[(index_i[e] * N + index_j[e]) + (sampleIndex * N * N)] > 0.01) {
				//printf("y[%d][%d]=%f\n", index_i[e], index_j[e], ySample[(index_i[e] * N + index_j[e]) + (sampleIndex * N * N)]);
				countY[index_i[e]][index_j[e]]++;
				countY[index_j[e]][index_i[e]]++;
			}
		}
		qsort((CostSort*)Sorted, counter, sizeof(Sorted[0]), Comparevalue);

	}

	
	for (i = 0; i < N - 1; i++) {
		for (j = i + 1; j < N; j++) {
			countY[i][j] = 0;
			countY[j][i] = 0;
		}
	}

	Out = open_file(outfile, "a+");
	fprintf(Out, "sampleIndex = %d\t", sampleIndex);
	fclose(Out);
	return 0;
}

int Comparevalue(const void* a, const void* b)
{
	if (((CostSort*)a)->value > ((CostSort*)b)->value)
		return 1;
	if (((CostSort*)a)->value < ((CostSort*)b)->value)
		return -1;
	return 0;
}
