#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

double **A, **L, **U;

void displayALU(int n) {
	printf("---------------A-----------------\n");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			printf("%0.12lf ", A[i][j]);
		}
		printf("\n");
	}
	printf("---------------L-----------------\n");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			printf("%0.12lf ", L[i][j]);
		}
		printf("\n");
	}
	printf("---------------U-----------------\n");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			printf("%0.12lf ", U[i][j]);
		}
		printf("\n");
	}
}

void displayLU(int n) {
	printf("---------------LU-----------------\n");
	double sum;
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			sum = 0.0;
			for (int k = 0; k < n; k++) {
				sum += L[i][k] * U[k][j];
			}
			printf("%0.12lf ", sum);
		}
		printf("\n");
	}
}

void strategy4(int n, int t) {

}

void get_filename(int n, int num_threads, int strategy, char *outfile, char *file) {
	char temp_thread[5], temp_strategy[5];
	sprintf(temp_thread, "%d", num_threads);
	sprintf(temp_strategy, "%d", strategy);
	strcpy(outfile, file);
	strcat(outfile, temp_strategy);
	strcat(outfile, "_");
	strcat(outfile, temp_thread);
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
		U[i][i]=1.0;
	}
}

int main(int argc, char *argv[]) {
	int n = atoi(argv[1]);
	char *input_file = argv[2];
	int comm_sz;
	int rank;

	init_matrix(n, input_file);

	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(rank==0){
		displayALU(n);
	}

	int i,k;
	double sum=0;
	int diff,start,end;

	for (int j = 0; j < n; j++) {
		diff = (n-j)/comm_sz;
		start = (diff)*rank+j; 
		end = (diff)*(rank+1)+j;
		if(rank==comm_sz-1){
			end = n;
		}
		for (i = start; i < end; i++) {
			sum = 0;
			for (k = 0; k < start; k++) {
				sum = sum + L[i][k] * U[k][j];	
			}
			L[i][j] = A[i][j] - sum;
		}

		for (i = start; i < end; i++) {
			sum = 0;
			for(k = 0; k < start; k++) {
				sum = sum + L[j][k] * U[k][i];
			}
			if (L[j][j] == 0) {	
				printf("Exiting!\n");			
			}
			U[j][i] = (A[j][i] - sum) / L[j][j];
		}		
		
	}



	MPI_Barrier(MPI_COMM_WORLD);
	
	if(rank==0){
		displayALU(n);
		displayLU(n);
		write_to_file(n, comm_sz, 4);
	}
	MPI_Finalize();
	return 0;
}