#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

double **A, **L, **U;

void matrix_mult(double **M1, double **M2, int n) {
	double sum;
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			sum = 0.0;
			for (int k = 0; k < n; k++) {
				sum += M1[i][k] * M2[k][j];
			}
			printf("%0.12lf ", sum);
		}
		printf("\n");
	}
}

void print_matrix(double **M, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			printf("%0.12lf ", M[i][j]);
		}
		printf("\n");
	}
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
}

void transpose_matrix(double **M, double **M_t, int n) {
	for(int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			//M_t[i*n + j] = M[j*n + i];
			M_t[i][j] = M[j][i];
		}
	}
}

void get_all_columns(double **M, double **col, int n, int rank, int cols_per_thread) {
	printf("\nGet All columns\n");
	for(int sub_col = 0; sub_col < cols_per_thread; sub_col++) {
		for(int row_index = 0; row_index < n; row_index++) {
			//col[sub_col*n + row_index] = M[(row_index * n) + (rank * cols_per_thread) + sub_col];
			col[sub_col][row_index] = M[row_index][rank*cols_per_thread + sub_col];
		}
	}
	printf("Got All columns\n\n");
}

void get_all_rows(double **M, double **row, int n, int rank, int rows_per_thread) {
	printf("\nGet All rows\n");
	for(int sub_row = 0; sub_row < rows_per_thread; sub_row++) {
		for(int col_index = 0; col_index < n; col_index++) {
			//row[sub_row*n + col_index] = M[(rank * rows_per_thread * n) + (sub_row * n) + col_index];
			row[sub_row][col_index] = M[rank*rows_per_thread + sub_row][col_index];
		}
	}
	printf("Got All rows\n\n");
}

void zero_fill(double **M, int n) {
	for(int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			//M[i*n + j] = 0;
			M[i][j] = 0;
		}
	}
}

void init_comp_matrix(double **M, int m , int n) {
	M = (double **)malloc(n * sizeof(double *));
	for (int i = 0; i < m; i++) {
		M[i] = (double *)malloc(n * sizeof(double));
	}
}

void copy_matrix(double **M_copy, double **M, int n) {
	for(int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			M_copy[i][j] = M[i][j];
		}
	}
} 

void LU_decomp(double **A, double **L, double **U, int n, int rank, int comm_sz) {
	int rows_per_thread = n/comm_sz;
	double **A_2, **L_2, **U_2, **L_3, **U_3, **U_t;

	A_2 = (double **)malloc(n * sizeof(double *));
	L_2 = (double **)malloc(n * sizeof(double *));
	U_2 = (double **)malloc(n * sizeof(double *));
	L_3 = (double **)malloc(n * sizeof(double *));
	U_3 = (double **)malloc(n * sizeof(double *));
	U_t = (double **)malloc(n * sizeof(double *));
	for (int i = 0; i < n; i++) {
		A_2[i] = (double *)malloc(n * sizeof(double));
		L_2[i] = (double *)malloc(n * sizeof(double));
		U_2[i] = (double *)malloc(n * sizeof(double));
		U_t[i] = (double *)malloc(n * sizeof(double));
	}
	for (int i = 0; i < rows_per_thread; i++) {
		L_3[i] = (double *)malloc(rows_per_thread * sizeof(double));
		U_3[i] = (double *)malloc(rows_per_thread * sizeof(double));
	}

	if (rank == 0) {
		copy_matrix(A_2, A, n);
		zero_fill(U_t, n);
	}
	zero_fill(L_2, n);
	zero_fill(U_2, n);
	MPI_Bcast(A_2, n*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	printf("ENTERING FOR LOOP!\n");

	for (int i = 0; i < n; i++) {
		printf("ENTERING FIRST INNER LOOP!\n");
		for (int j = rank*rows_per_thread; j < (rank + 1)*rows_per_thread; j++) {
			if (j < i) {
				//L_2[j*n + i] = 0;
				L_2[j][i] = 0;
			}
			else {
				L_2[j*n + i] = A_2[j*n + i];
				for (int k = 0; k < i; k++) {
					//L_2[j*n + i] = L_2[j*n + i] - L_2[j*n + k] * U_2[k*n + i];
					L_2[j][i] = L_2[j][i] - L_2[j][k] * U_2[k][i];
				}
			}
			
		}
		printf("LEAVING FIRST INNER LOOP!\n");
		get_all_rows(L_2, L_3, n, rank, rows_per_thread);
		MPI_Gather(L_3, n*rows_per_thread, MPI_DOUBLE, L_2, n*rows_per_thread, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(L_2, n*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		printf("ENTERING SECOND INNER LOOP!\n");
		for(int j = rank*rows_per_thread; j < (rank + 1)*rows_per_thread; j++) {
			if(j < i) {
				//U_2[i*n + i] = 0;
				U_2[i][i] = 0;
			}
			else if(j == i) {
				//U_2[i*n + j] = 1;
				U_2[i][j] = 1;
			}
			else {
				//U_2[i*n + j] = A_2[i*n + j] / L_2[i*n + i];
				U_2[i][j] = A_2[i][j] / L_2[i][i];
				for(int k = 0; k < i; k++) {
					//U_2[i*n + j] = U_2[i*n + j] - ((L_2[i*n + k] * U_2[k*n + j]) / L_2[i*n + i]);
					U_2[i][j] = U_2[i][j] - ((L_2[i][k] * U_2[k][j]) / L_2[i][i]);
				}
			}
		}
		printf("LEAVING SECOND INNER LOOP!\n");
		get_all_columns(U_2, U_3, n, rank, rows_per_thread);
		MPI_Gather(U_3, n*rows_per_thread, MPI_DOUBLE, U_t, n*rows_per_thread, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		transpose_matrix(U_t, U_2, n);
		MPI_Bcast(U_2, n*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	printf("LEAVING FOR LOOP!\n");
	if(rank == 0) {
		printf("---------------L-----------------\n");
		print_matrix(L_2, n);
		printf("---------------U-----------------\n");
		print_matrix(U_2, n);
		printf("---------------LU-----------------\n");
		matrix_mult(L_2, U_2, n);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

void ini(double M[], int n) {
	for(int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			M[n*i + j] = 2.0;
		}
	}
}

int main(int argc, char *argv[]) {
	int n = atoi(argv[1]);
	char *input_file = argv[2];

	double K[n*n];
	ini(K, n);
	for(int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			printf("%f ", K[i*n + j]);
		}
		printf("\n");
	}


	double time_taken;
	int comm_sz; // Number of Processes
	int rank; // Process number

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0) {
		init_matrix(n, input_file);
		printf("---------------A-----------------\n");
		print_matrix(A, n);
		time_taken = MPI_Wtime();
	}
	LU_decomp(A, L, U, n, rank, comm_sz);
	if(rank == 0) {
		time_taken = MPI_Wtime() - time_taken;
		printf("Time taken in strategy 4 is %f seconds\n", time_taken);
	}
	MPI_Finalize();
	return 0;
}