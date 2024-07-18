#include "headers.h"

int** interdiction, ** fortification, ** auxiliaryInterdiction, ** incumbantInterdiction, ** incumbantfortification, ** incumbantY, ** yAuxiliary, ** interdictionperceived, ** InterdictionS, ** waitingfortification, ** waitingInterdiction;
int Fbudget, *Ibudget, Fortification_Flag, CoverConstraint, FlagCover, initialsample, S;
int sample, sampleIndex, bindex, iteration, iter, sampleRemoved, acceptedSample, sample_index;
double* objSample, * ySample, zPerceived, zbar, M, upperboundsto, lowerboundsto, * p, Ubound, ratio, *waitinglistRealdamage;
int checkVisted, tempN, sglobal, globalStop, waitinglist, indexGlobal, flag;
double incumbantzbar, incumbantUB, incumbantLB;
int** countY, MaxCountY;
char OutName[100];
FILE* Out = NULL;
char outfile[20];
double epsilonBSF;

double Timelimit = 24 * 3600;
double SampleTimeLimit = 10 * 3600;
double checktime, TotalTime, SampleTime, timeIncumbent;

////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) {

	char	instance[20];
	char	path[50];
	FILE* ini;
	clock_t  start, end;

	GLOBAL_INFO global;
	CONFIGURATION param;		// Global Paramenters of the algorithm
	INSTANCE* instances;		// Array of instance data
	RESULTS* results;			// Array of results
	int n_of_instances;			// Number of instances to solve
	int i, e, ii, bindex, j, index1, index;						// Index variable

	int status = 0;

	if (argc < 3) { show_help(argv[0]); }
	else {
		// read experimental settings:
		status = read_CONFIGURATION(argv[1], &param);
		assert(status == 0);

		// read the instance list
		instances = create_INSTANCE_vector(&n_of_instances, argv[2]);
		results = create_RESULTS_vector(n_of_instances);

		timeNow = time(NULL);
		tm = localtime(&timeNow);
		Out = open_file(outfile, "a+");
		fprintf(Out, "\n\n %s\n", asctime(tm));
		fclose(Out);

		for (i = 0; i < n_of_instances; i++) {
			double TimeBSF, timewhile, timenew;
			//double startSampling, startTotal, TimeBSF, endBSF, timewhile, timenew;//
			//double endSampling, endTotal, endincumbent;//
			TotalTime = 0;
			//startSampling = 0;//
			SampleTime = 0;
			//startTotal = 0;//
			TotalTime = 0;
			//SampleTime = 0;//
			//endSampling = 0;//
			//endTotal = 0;//
			//endBSF = 0;//
			//double startwhile, endwhile, startnew, endnew;//

			countfails = 0;
			// Read the instance data:
			status = read_INSTANCE(&(instances[i]));
			assert(status == 0);

			int counterImprovement = 0;
			double zbarPrevious;

			iter = 0;
			
			zbar = CPX_INFBOUND;

			//startTotal = clock();//
			//startSampling = clock();//
			gettimeofday(&startTotal, NULL);
			gettimeofday(&startSampling, NULL);

			status = sampling(param, instances[i], &(results[i]));

			gettimeofday(&endSampling, NULL);
			SampleTime = ((double)(endSampling.tv_sec - startSampling.tv_sec) * 1000 + (double)(endSampling.tv_usec - startSampling.tv_usec) / 1000) / 1000;


			//endSampling = clock();//
			//SampleTime = (double)(endSampling - startSampling) / (double)(CLOCKS_PER_SEC);//

			TimeBSF = 0;
			CoverConstraint = 0;

			for (ii = 0; ii < tempN - 1; ii++) {
				for (j = ii + 1; j < tempN; j++) {
					fortification[ii][j] = 0;
				}
			}
			int iter = 0;
			
			incumbantzbar = 0;
			incumbantUB = 0;
			incumbantLB = 0;
			Fortification_Flag = 0;
			
			globalStop = 0;


			assert(status == 0);
			while (TimeBSF < Timelimit && globalStop == 0 && sampleIndex < iteration) {
				while (Fortification_Flag == 0 && TimeBSF < Timelimit && sampleIndex < iteration) {
					FortificationProblem(param, instances[i], &(results[i]));
					if (Fortification_Flag == 1) {
						break;
					}
					upperboundsto = CPX_INFBOUND;
					lowerboundsto = 0;
					while (lowerboundsto < zbar && sampleIndex < iteration) {
						FlagCover = 0;

						zbarPrevious = zbar;

						upperboundsto = 0;
						lowerboundsto = 0;

						for (sglobal = 0; sglobal < S; sglobal++) {
							perceivedDamage(param, instances[i], &(results[i]));
							upperboundsto += p[sglobal] * Ubound;
							lowerboundsto += p[sglobal] * realDamage(param, instances[i], &(results[i]));

							sampleIndex++;

						}

						if (upperboundsto < zbar) {
							zbar = upperboundsto;
							for (sample_index = 0; sample_index < sampleIndex; sample_index++) {
								if (objSample[sample_index] > upperboundsto + 0.000001) { ///if obj>UB, we remove that sample
									for (index = sample_index; index < sampleIndex; index++) {
										objSample[index] = objSample[index + 1];
										for (ii = 0; ii < tempN - 1; ii++) {
											for (j = ii + 1; j < tempN; j++) {
												ySample[(ii * tempN + j) + (index * tempN * tempN)] = ySample[(ii * tempN + j) + ((index + 1) * tempN * tempN)];
											}
										}
									}
									sampleIndex--;
									sample_index--;
								}
							}
						}

						if (zbar - zbarPrevious == 0) {
							counterImprovement++;
						}

						if (lowerboundsto >= zbar) {
							PreFortification(param, instances[i], &(results[i]));
						}
						if (ABS(lowerboundsto - upperboundsto) < 0.000001 && ABS(upperboundsto - zbar) < 0.000001) {
							status = backwardAlgorithm(param, instances[i], &(results[i]));
							gettimeofday(&endincumbent, NULL);
							timeIncumbent = ((double)(endincumbent.tv_sec - startTotal.tv_sec) * 1000 + (double)(endincumbent.tv_usec - startTotal.tv_usec) / 1000) / 1000;

							//endincumbent = clock();//
							//timeIncumbent = (double)(endincumbent - startTotal) / (double)(CLOCKS_PER_SEC);//

						}

						if (counterImprovement > 30) {
							counterImprovement = 0;
							storeWaitList(param, instances[i], &(results[i]));
							lowerboundsto = CPX_INFBOUND;
						}

						ratio = (zbar - lowerboundsto) / zbar;
						if (ratio <= epsilonBSF && lowerboundsto < zbar) {
							storeWaitList(param, instances[i], &(results[i]));
							lowerboundsto = CPX_INFBOUND;
						}

						gettimeofday(&endBSF, NULL);
						TimeBSF = ((double)(endBSF.tv_sec - startTotal.tv_sec) * 1000 + (double)(endBSF.tv_usec - startTotal.tv_usec) / 1000) / 1000;

						//endBSF = clock(); //
						//TimeBSF = (double)(endBSF - startTotal) / (double)(CLOCKS_PER_SEC);//

						printf("LB = %.3f \t UB =  %.3f \t Zbar =  %.3f/n", lowerboundsto, upperboundsto, zbar);
					}
				}

				if (waitinglist == 0) {
					globalStop = 1;
				}

				else if (TimeBSF < Timelimit && sampleIndex < iteration) {
					Fortification_Flag = 0;
					epsilonBSF = 0;
					for (indexGlobal = 0; indexGlobal < waitinglist; indexGlobal++) {
						if (waitinglistRealdamage[indexGlobal] > zbar) {
							addCoverconstraint(param, instances[i], &(results[i]));
							for (index = indexGlobal; index < waitinglist; index++) {
								waitinglistRealdamage[index] = waitinglistRealdamage[index + 1];
								for (sglobal = 0; sglobal < S; sglobal++) {
									for (ii = 0; ii < tempN - 1; ii++) {
										for (j = ii + 1; j < tempN; j++) {
											waitingInterdiction[(ii * tempN + j) + (sglobal * tempN * tempN)][index] = waitingInterdiction[(ii * tempN + j) + (sglobal * tempN * tempN)][index + 1];
										}
									}
								}
							}
							waitinglist--;
							indexGlobal--;
						}
					}
					for (indexGlobal = 0; indexGlobal < waitinglist; indexGlobal++) {
						flag = checkfeasibility();
						if (flag == 1) {
							for (sglobal = 0; sglobal < S; sglobal++) {
								for (ii = 0; ii < tempN - 1; ii++) {
									for (j = ii + 1; j < tempN; j++) {
										fortification[ii][j] = waitingfortification[(ii * tempN + j) + (sglobal * tempN * tempN)][indexGlobal];
									}
								}
							}
							lowerboundsto = waitinglistRealdamage[indexGlobal];
							while (lowerboundsto < zbar && TimeBSF < Timelimit && sampleIndex < iteration) {
								FlagCover = 0;
								upperboundsto = 0;
								lowerboundsto = 0;

								for (sglobal = 0; sglobal < S; sglobal++) {
									perceivedDamage(param, instances[i], &(results[i]));
									upperboundsto += p[sglobal] * Ubound;
									lowerboundsto += p[sglobal] * realDamage(param, instances[i], &(results[i]));
									sampleIndex++;
								}
								if (upperboundsto < zbar) {
									zbar = upperboundsto;
									timenew = 0;
									for (sample_index = 0; sample_index < sampleIndex; sample_index++) {
										if (objSample[sample_index] > upperboundsto + 0.000001) { ///if obj>UB, we remove that sample
											for (index = sample_index; index < sampleIndex; index++) {
												objSample[index] = objSample[index + 1];
												for (ii = 0; ii < tempN - 1; ii++) {
													for (j = ii + 1; j < tempN; j++) {
														ySample[(ii * tempN + j) + (index * tempN * tempN)] = ySample[(ii * tempN + j) + ((index + 1) * tempN * tempN)];
													}
												}
											}
											sampleIndex--;
											sample_index--;
										}
									}
								}

								if (lowerboundsto >= zbar) {
									PreFortification(param, instances[i], &(results[i]));
								}
								if (ABS(lowerboundsto - upperboundsto) < 0.000001 && ABS(upperboundsto - zbar) < 0.000001) {
									status = backwardAlgorithm(param, instances[i], &(results[i]));
									gettimeofday(&endincumbent, NULL);
									timeIncumbent = ((double)(endincumbent.tv_sec - startTotal.tv_sec) * 1000 + (double)(endincumbent.tv_usec - startTotal.tv_usec) / 1000) / 1000;

								}
								ratio = (zbar - lowerboundsto) / zbar;
								if (ratio <= epsilonBSF && lowerboundsto < zbar) {
									storeWaitList(param, instances[i], &(results[i]));
									lowerboundsto = CPX_INFBOUND;
								}
								gettimeofday(&endBSF, NULL);
								TimeBSF = ((double)(endBSF.tv_sec - startTotal.tv_sec) * 1000 + (double)(endBSF.tv_usec - startTotal.tv_usec) / 1000) / 1000;

								//endBSF = clock();//
								//TimeBSF = (double)(endBSF - startTotal) / (double)(CLOCKS_PER_SEC);//
							}
						}

						for (index = indexGlobal; index < waitinglist; index++) {
							waitinglistRealdamage[index] = waitinglistRealdamage[index + 1];
							for (sglobal = 0; sglobal < S; sglobal++) {
								for (ii = 0; ii < tempN - 1; ii++) {
									for (j = ii + 1; j < tempN; j++) {
										waitingInterdiction[(ii * tempN + j) + (sglobal * tempN * tempN)][index] = waitingInterdiction[(ii * tempN + j) + (sglobal * tempN * tempN)][index + 1];
									}
								}
							}

						}
						waitinglist--;
						indexGlobal--;
					}
					epsilonBSF = 0.25;
				}



			}
			////////////////////////////////////////////////////////////////////
			status = backwardAlgorithmfinal(param, instances[i], &(results[i]));
		
			gettimeofday(&endTotal, NULL);
			TotalTime = ((double)(endTotal.tv_sec - startTotal.tv_sec) * 1000 + (double)(endTotal.tv_usec - startTotal.tv_usec) / 1000) / 1000;

			//endTotal = clock(); //
			//TotalTime = (double)(endTotal - startTotal) / (double)(CLOCKS_PER_SEC);//

			printf("/nFinished 2");
			Out = open_file(outfile, "a+");
			fprintf(Out, "Sample time = %f \t", SampleTime);
			fprintf(Out, "Time incumbent = %f \t", timeIncumbent);
			fprintf(Out, "Total time = %f \t", TotalTime);
			fprintf(Out, "\n");

			fclose(Out);


			// Free instance data:
			free_INSTANCE(&(instances[i]));
			////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////
		}

		// Output average times & gaps...
		//output_averaged_RESULTS(filename, param, n_of_instances, results);

		// clean:
		//free_RESULTS_vector(&results, n_of_instances);
		free_INSTANCE_vector(&instances);
	}

	return status; // 0 meaning success
}
////////////////////////////////////////////////////////////////////////////////
int checkfeasibility() {
	int ii, j, k;
	for (ii = 0; ii < CoverConstraint; ii++) {
		int sum = 0;
		for (sglobal = 0; sglobal < S; sglobal++) {
			for (k = 0; k < tempN - 1; k++) {
				for (j = k + 1; j < tempN; j++) {
					sum += auxiliaryInterdiction[k * tempN + j][ii * S + sglobal] * waitingfortification[(k * tempN + j) + (sglobal * tempN * tempN)][indexGlobal];
				}
			}
		}
		/*if (sum  == Ibudget) {
			flag == 1;
		}*/
		if (sum < 1) {
			return 0;
		}
	}
	return 1;

}