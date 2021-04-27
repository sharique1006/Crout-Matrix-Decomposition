#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

double **A, **L, **U;

void sequential(int n) {
	int i, j, k;
	double sum = 0;
	L = (double **)malloc(n * sizeof(double *));
	U = (double **)malloc(n * sizeof(double *));
	for(int i = 0; i < n; i++) {
		L[i] = (double *)malloc(n * sizeof(double *));
		U[i] = (double *)malloc(n * sizeof(double *));
		for (int j = 0; j < n; j++) {
			L[i][j] = 0.0;
			U[i][j] = 0.0;
		}
		U[i][i] = 1;
	}

	// for (i = 0; i < n; i++) {
	// 	U[i][i] = 1;
	// }

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
	L = (double **)malloc(n * sizeof(double *));
	U = (double **)malloc(n * sizeof(double *));
	#pragma omp parallel for
	for(int i = 0; i < n; i++) {
		L[i] = (double *)malloc(n * sizeof(double *));
		U[i] = (double *)malloc(n * sizeof(double *));
		for (int j = 0; j < n; j++) {
			L[i][j] = 0.0;
			U[i][j] = 0.0;
		}
		U[i][i] = 1;
	}
	// double sum = 0;
	// omp_set_nested(1);

	// #pragma omp parallel for
	// for (int i = 0; i < n; i++) {
	// 	U[i][i] = 1;
	// }
	for (int j = 0; j < n; j++) {
		#pragma omp parallel for schedule(static)
		for (int i = j; i < n; i++) {
			double sum = 0;
			for (int k = 0; k < j; k++) {
				sum = sum + L[i][k] * U[k][j];	
			}
			L[i][j] = A[i][j] - sum;
		}
		if (L[j][j] == 0) {	
			printf("Exiting!\n");			
			exit(0);
		}
		#pragma omp parallel for schedule(static)
		for (int i = j; i < n; i++) {
			double sum = 0;
			for(int k = 0; k < j; k++) {
				sum = sum + L[j][k] * U[k][i];
			}
			
			U[j][i] = (A[j][i] - sum) / L[j][j];
		}
	}
}

void strategy2(int n, int t) {
	omp_set_num_threads(t);

	// // #pragma omp parallel for
	// for (int i = 0; i < n; i++) {
	// 	U[i][i] = 1;
	// }
	L = (double **)malloc(n * sizeof(double *));
	U = (double **)malloc(n * sizeof(double *));
	for(int i = 0; i < n; i++) {
		L[i] = (double *)malloc(n * sizeof(double *));
		U[i] = (double *)malloc(n * sizeof(double *));
		for (int j = 0; j < n; j++) {
			L[i][j] = 0.0;
			U[i][j] = 0.0;
		}
		U[i][i] = 1;
	}
	// #pragma omp parallel
	for (int j = 0; j < n; j++) {
		double sum = 0;
		for (int k = 0; k < j; k++) {
			sum = sum + L[j][k] * U[k][j];	
		}
		L[j][j] = A[j][j] - sum;
		if (L[j][j] == 0) {	
			printf("Exiting!\n");			
			exit(0);
		}
		#pragma omp parallel sections
		{
			#pragma omp section
			{
				#pragma omp parallel for schedule(static)
				for (int i = j+1; i < n; i++) {
					double sum = 0;
					for (int k = 0; k < j; k++) {
						sum = sum + L[i][k] * U[k][j];	
					}
					L[i][j] = A[i][j] - sum;
				}
			}

			#pragma omp section
			{	
				#pragma omp parallel for schedule(static)
				for (int i = j; i < n; i++) {
					double sum = 0;
					for(int k = 0; k < j; k++) {
						sum = sum + L[j][k] * U[k][i];
					}
					
					U[j][i] = (A[j][i] - sum) / L[j][j];
				}
			}
		}
	}
}

void strategy3(int n, int num_threads) {
	omp_set_num_threads(num_threads);
	L = (double **)malloc(n * sizeof(double *));
	U = (double **)malloc(n * sizeof(double *));
	#pragma omp parallel for
	for(int i = 0; i < n; i++) {
		L[i] = (double *)malloc(n * sizeof(double *));
		U[i] = (double *)malloc(n * sizeof(double *));
		for (int j = 0; j < n; j++) {
			L[i][j] = 0.0;
			U[i][j] = 0.0;
		}
		U[i][i] = 1;
	}
	for (int j = 0; j < n; j++) {
		double sum = 0;
		// #pragma omp parallel for schedule(dynamic)
		for (int k = 0; k < j; k++) {
			sum = sum + L[j][k] * U[k][j];	
		}
		L[j][j] = A[j][j] - sum;
		if (L[j][j] == 0) {	
			printf("Exiting!\n");			
			exit(0);
		}
		#pragma omp parallel sections
		{
			#pragma omp section
			{
				for (int i = j+1; i < n; i++) {
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
					
					U[j][i] = (A[j][i] - sum) / L[j][j];
				}
			}
		}
	}
}

void write_output(char *fname, double** arr, int n ){
	FILE *f = fopen(fname, "w");
	for( int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			fprintf(f, "%0.12f ", arr[i][j]);
		}
		fprintf(f, "\n");
	}
	fclose(f);
}

void init_matrix(int n, char *input_file) {
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
}

int main(int argc, char *argv[]) {
	int n = atoi(argv[1]);
	char *input_file = argv[2];
	int num_threads = atoi(argv[3]);
	int strategy = atoi(argv[4]);
	char *outfile_L = argv[5];
	char *outfile_U = argv[6];

	init_matrix(n, input_file);

	if (strategy == 0) sequential(n);
	else if (strategy == 1) strategy1(n, num_threads);
	else if (strategy == 2) strategy2(n, num_threads);
	else if (strategy == 3) strategy3(n, num_threads);

	write_output(outfile_L, L, n);
	write_output(outfile_U, U, n);

	return 0;
}