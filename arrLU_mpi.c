#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

void get_filename(int num_threads, int strategy, char *outfile, char *file) {
	char temp_thread[5], temp_strategy[5];
	sprintf(temp_thread, "%d", num_threads);
	sprintf(temp_strategy, "%d", strategy);
	strcpy(outfile, file);
	strcat(outfile, temp_strategy);
	strcat(outfile, "_");
	strcat(outfile, temp_thread);
	strcat(outfile, ".txt");
}

void write_output(double arr_L[],  double arr_U[], int n, int num_threads, int strategy) {
	char outfile_L[100], outfile_U[100];
	get_filename(num_threads, strategy, outfile_L, "output_L_");
	get_filename(num_threads, strategy, outfile_U, "output_U_");
	FILE *f1 = fopen(outfile_L, "w");
	FILE *f2 = fopen(outfile_U, "w");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			fprintf(f1, "%0.12lf ", arr_L[n*i + j]);
			fprintf(f2, "%0.12lf ", arr_U[n*i + j]);
		}
		fprintf(f1, "\n");
		fprintf(f2, "\n");
	}
	fclose(f1);
	fclose(f2);
}


void init_matrix(double A[], int n, char *input_file) {
	FILE *f;
	f = fopen(input_file, "r");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (!fscanf(f, "%140lf ", &A[n*i + j])) {
				break;
			}
		}
	}
	fclose(f);
}

void print_matrix(double M[], int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			printf("%0.12lf ", M[n*i + j]);
		}
		printf("\n");
	}
}

void matrix_mult(double M1[], double M2[], int n) {
	double sum;
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			sum = 0.0;
			for (int k = 0; k < n; k++) {
				sum += M1[n*i + k] * M2[n*k + j];
			}
			printf("%0.12lf ", sum);
		}
		printf("\n");
	}
}

void transpose_matrix(double M[], double M_t[], int n) {
	for(int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			M_t[i*n + j] = M[j*n + i];
		}
	}
}

void get_all_columns(double M[], double col[], int n, int rank, int cols_per_thread) {
	for(int sub_col = 0; sub_col < cols_per_thread; sub_col++) {
		for(int row_index = 0; row_index < n; row_index++) {
			col[sub_col*n + row_index] = M[(row_index * n) + (rank * cols_per_thread) + sub_col];
		}
	}
}

void get_all_rows(double M[], double row[], int n, int rank, int rows_per_thread) {
	for(int sub_row = 0; sub_row < rows_per_thread; sub_row++) {
		for(int col_index = 0; col_index < n; col_index++) {
			row[sub_row*n + col_index] = M[(rank * rows_per_thread * n) + (sub_row * n) + col_index];
		}
	}
}

void zero_fill(double M[], int n) {
	for(int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			M[i*n + j] = 0;
		}
	}
}

void copy_matrix(double M_copy[], double M[], int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			M_copy[n*i + j] = M[n*i + j];
		}
	}
}

void LU_decomp(double A[], double L[], double U[], int n, int rank, int comm_sz) {
	int rows_per_thread = n/comm_sz;

	double A_2[n*n], L_2[n*n], U_2[n*n], L_3[n*rows_per_thread], U_3[n*rows_per_thread],  U_t[n*n];

	if (rank == 0) {
		copy_matrix(A_2, A, n);
		zero_fill(U_t, n);
	}
	zero_fill(L_2, n);
	zero_fill(U_2, n);
	MPI_Bcast(A_2, n*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// printf("ENTERING FOR LOOP!\n");

	for (int i = 0; i < n; i++) {
		// printf("ENTERING FIRST INNER LOOP!\n");
		for (int j = rank*rows_per_thread; j < (rank + 1)*rows_per_thread; j++) {
			if (j < i) {
				L_2[j*n + i] = 0;
			}
			else {
				L_2[j*n + i] = A_2[j*n + i];
				for (int k = 0; k < i; k++) {
					L_2[j*n + i] = L_2[j*n + i] - L_2[j*n + k] * U_2[k*n + i];
				}
			}
			
		}
		// printf("LEAVING FIRST INNER LOOP!\n");
		get_all_rows(L_2, L_3, n, rank, rows_per_thread);
		MPI_Gather(L_3, n*rows_per_thread, MPI_DOUBLE, L_2, n*rows_per_thread, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(L_2, n*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		// printf("ENTERING SECOND INNER LOOP!\n");
		for(int j = rank*rows_per_thread; j < (rank + 1)*rows_per_thread; j++) {
			if(j < i) {
				U_2[i*n + i] = 0;
			}
			else if(j == i) {
				U_2[i*n + j] = 1;
			}
			else {
				U_2[i*n + j] = A_2[i*n + j] / L_2[i*n + i];
				for(int k = 0; k < i; k++) {
					U_2[i*n + j] = U_2[i*n + j] - ((L_2[i*n + k] * U_2[k*n + j]) / L_2[i*n + i]);
				}
			}
		}
		// printf("LEAVING SECOND INNER LOOP!\n");
		get_all_columns(U_2, U_3, n, rank, rows_per_thread);
		MPI_Gather(U_3, n*rows_per_thread, MPI_DOUBLE, U_t, n*rows_per_thread, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		transpose_matrix(U_t, U_2, n);
		MPI_Bcast(U_2, n*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	// printf("LEAVING FOR LOOP!\n");
	if(rank == 0) {
		printf("---------------L-----------------\n");
		print_matrix(L_2, n);
		printf("---------------U-----------------\n");
		print_matrix(U_2, n);
		printf("---------------LU-----------------\n");
		matrix_mult(L_2, U_2, n);
		write_output(L_2, U_2, n, comm_sz, 4);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

int main(int argc, char *argv[]) {
	int n = atoi(argv[1]);
	char *input_file = argv[2];

	double A[n*n], L[n*n], U[n*n];

	double time_taken;
	int comm_sz; // Number of Processes
	int rank; // Process number

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0) {
		init_matrix(A, n, input_file);
		printf("---------------A-----------------\n");
		print_matrix(A, n);
		time_taken = MPI_Wtime();
	}
	LU_decomp(A, L, U, n, rank, comm_sz);
	if(rank == 0) {
		time_taken = MPI_Wtime() - time_taken;
		printf("\nTime taken in strategy 4 is %f seconds\n", time_taken);
	}
	MPI_Finalize();
	return 0;
}