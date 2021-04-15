#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

double **A, **L, **U;

void sequential(int n) {
	int i, j, k;
	double sum = 0;

	for (i = 0; i < n; i++) {
		U[i][i] = 1;
	}

	for (j = 0; j < n; j++) {
		for (i = j; i < n; i++) {
			sum = 0;
			for (k = 0; k < j; k++) {
				sum = sum + L[i][k] * U[k][j];	
			}
			L[i][j] = A[i][j] - sum;
		}

		for (i = j; i < n; i++) {
			sum = 0;
			for(k = 0; k < j; k++) {
				sum = sum + L[j][k] * U[k][i];
			}
			if (L[j][j] == 0) {	
				printf("Exiting!\n");			
				exit(0);
			}
			U[j][i] = (A[j][i] - sum) / L[j][j];
		}
	}
}

void strategy1(int n, int num_threads) {
	omp_set_num_threads(num_threads);
	double sum = 0;
	// omp_set_nested(1);

	#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		U[i][i] = 1;
	}
	// #pragma omp parallel shared(A,L,U)
		for (int j = 0; j < n; j++) {
			#pragma omp parallel for schedule(dynamic)
			for (int i = j; i < n; i++) {
				sum = 0;
				L[i][j] = A[i][j];
				for (int k = 0; k < j; k++) {
					// sum = sum + L[i][k] * U[k][j];	
					L[i][j] = L[i][j] - L[i][k] * U[k][j];
				}
				// L[i][j] = A[i][j] - sum;
			}
			#pragma omp parallel for schedule(dynamic)
			for (int i = j; i < n; i++) {
				sum = 0;
				U[j][i] = A[j][i] / L[j][j];
				for(int k = 0; k < j; k++) {
					U[j][i] = U[j][i] - ((L[j][k] * U[k][i]) / L[j][j]);
				}
				if (L[j][j] == 0) {	
					printf("Exiting!\n");			
					exit(0);
				}
				// U[j][i] = (A[j][i] - sum) / L[j][j];
			}
		}
}

void strategy2(int n, int t) {
	for (int i = 0; i < n; i++) {
		U[i][i] = 1;
	}
	omp_set_num_threads(t);

	#pragma omp parallel
	for (int j = 0; j < n; j++) {
		#pragma omp parallel sections
		{
			#pragma omp section
			{
				for (int i = j; i < n; i++) {
					double sum = 0;
					for (int k = 0; k < j; k++) {
						sum = sum + L[i][k] * U[k][j];	
					}
					L[i][j] = A[i][j] - sum;
				}
			}
			#pragma omp section
			{
				for (int i = j; i < n; i++) {
					double sum = 0;
					for(int k = 0; k < j; k++) {
						sum = sum + L[j][k] * U[k][i];
					}
					if (L[j][j] == 0) {	
						printf("Exiting!\n");			
						exit(0);
					}
					U[j][i] = (A[j][i] - sum) / L[j][j];
				}
			}
		}
	}
}

void strategy3(int n, int num_threads) {

}

void get_filename(int n, int num_threads, int strategy, char *outfile, char *file) {
	char temp_thread[5], temp_strategy[5];
	sprintf(temp_thread, "%d", num_threads);
	sprintf(temp_strategy, "%d", strategy);
	strcpy(outfile, file);
	strcat(outfile, temp_thread);
	strcat(outfile, "_");
	strcat(outfile, temp_strategy);
	strcat(outfile, ".txt");
}

void write_to_file(int n, int num_threads, int strategy) {
	char outfile_L[100], outfile_U[100];
	get_filename(n, num_threads, strategy, outfile_L, "output_L_");
	get_filename(n, num_threads, strategy, outfile_U, "output_U_");
	FILE *f1 = fopen(outfile_L, "w");
	FILE *f2 = fopen(outfile_U, "w");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			fprintf(f1, "%0.12lf ", L[i][j]);
			fprintf(f2, "%0.12lf ", U[i][j]);
		}
		fprintf(f1, "\n");
		fprintf(f2, "\n");
	}
	fclose(f1);
	fclose(f2);
}

void init_matrix(int n, char *input_file) {
	/*Read A from file*/
	FILE *f;
	f = fopen(input_file, "r");
	A = (double **)malloc(n * sizeof(double *));
	for (int i = 0; i < n; i++) {
		A[i] = (double *)malloc(n * sizeof(double));
		for (int j = 0; j < n; j++) {
			if (!fscanf(f, "%140lf ", &A[i][j])) {
				break;
			}
		}
	}
	fclose(f);
	/*Initialise L & U as 0's*/
	L = (double **)malloc(n * sizeof(double *));
	U = (double **)malloc(n * sizeof(double *));
	for(int i = 0; i < n; i++) {
		L[i] = (double *)malloc(n * sizeof(double *));
		U[i] = (double *)malloc(n * sizeof(double *));
		for (int j = 0; j < n; j++) {
			L[i][j] = 0.0;
			U[i][j] = 0.0;
		}
	}
}

int main(int argc, char *argv[]) {
	int n = atoi(argv[1]);
	char *input_file = argv[2];
	int num_threads = atoi(argv[3]);
	int strategy = atoi(argv[4]);

	init_matrix(n, input_file);

	double start, end;

	if (strategy == 0) { start = omp_get_wtime(); sequential(n); end = omp_get_wtime(); }
	else if (strategy == 1) { start = omp_get_wtime(); strategy1(n, num_threads); end = omp_get_wtime(); }
	else if (strategy == 2) { start = omp_get_wtime(); strategy2(n, num_threads); end = omp_get_wtime(); }
	else if (strategy == 3) { start = omp_get_wtime(); strategy3(n, num_threads); end = omp_get_wtime(); }

	double time_taken = (end - start);
	printf("Time taken in strategy %d is %f seconds\n", strategy, time_taken);
	write_to_file(n, num_threads, strategy);

	return 0;
}