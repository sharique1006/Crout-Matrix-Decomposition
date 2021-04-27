#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

double *A;
float diff;
int j, n, comm_sz;

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
	// sprintf(outfile_L,"output_L.txt");
	// sprintf(outfile_U,"output_U.txt");

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

void init_matrix(int n, char *input_file) {
	
	FILE *f;
	f = fopen(input_file, "r");
	A = (double *)malloc(n*n * sizeof(double ));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			// if (!fscanf(f, "%140lf ", &A[i*n+j])) {
			if (!fscanf(f, "%lf ", &A[i*n+j])) {
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

void zero_fill(double *M, int n) {
	for(int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			M[i*n + j] = 0;
		}
	}
}

void zero_fill2(double *M, int n) {
	for(int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			M[i*n + j] = 0;
		}
		M[i*n + i] = 1;
	}
}

void copy_matrix(double *M_copy, double *M, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			M_copy[n*i + j] = M[n*i + j];
		}
	}
}

void Recv(void* location, int size, MPI_Datatype datatype, int source, int tag){
	if(size==-1){
		int rstart = diff * source + j;
		// int rend = diff * (source + 1) + j;
		int rend = (source==comm_sz-1)? n : diff * (source+1) +j;
		size = rend - rstart;
	}

	int tot = 0;
	int count;
	while(tot!=size){
		MPI_Status status;
		MPI_Recv(location,size, datatype, source, tag, MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, MPI_DOUBLE, &count);
		tot+=count;
	}
}

void LU_decomp(int n, int rank, int comm_sz) {

	double *L=(double *)malloc(n*n * sizeof(double));
	double *U_2=(double *)malloc(n*n * sizeof(double));

	for(int i=0;i<n;i++){
		for(int k=0;k<n;k++){
			L[i*n+k]=0;
			U_2[i*n+k]=0;
		}
	}

	for(int i = 0; i < n; i++) {
		U_2[i*n + i] = 1;
	}
	int i,k;
	int start,end;

	// if(rank==0){
	// 	for(int proc = 1; proc< comm_sz; proc++) MPI_Send(U_2,n*n,MPI_DOUBLE,proc,1,MPI_COMM_WORLD);
	// }else{
	// 	Recv(U_2, n*n, MPI_DOUBLE, 0, 1);
	// }

	for (j = 0; j < n; j++) {
		diff = (float)(n-j)/(float)(comm_sz);
		start = (diff)*rank+j; 
		end = (rank==comm_sz-1)? n : (diff)*(rank+1)+j;

		// Calc L
		if(rank==0){
			// printf("=====================================\n");
			// printf("Thread %d Receving %d size\n",rank,end-start);
			for (i = start; i < end; i++) {
				double currSum = 0;
				for (k = 0; k < j; k++) {
					currSum = currSum + L[i*n+k] * U_2[k*n+j];
					// printf("Mult L[%d,%d](%lf) and U[%d,%d](%lf) for L[%d,%d] with A[%d][%d](%lf)\n",i,k,L[i*n+k],k,j,U_2[k*n+j],i,j,i,j,A[(i*n+j)]);
				}
				L[i*n+j] = A[i*n+j] - currSum;
				// printf("Thread %d wrote in L[%d,%d]\n",rank,i,j);
			}

			for(int proc = 1; proc<comm_sz; proc++){
				int rstart = diff * proc + j;
				int rend = (proc==comm_sz-1)? n : diff * (proc+1) +j;
				double recsum[rend-rstart];
				Recv(recsum, -1, MPI_DOUBLE, proc, 0);

				for (i = rstart; i < rend; i++) {
					L[i*n+j] = A[i*n+j] - recsum[i-rstart];
					// printf("Thread %d wrote in L[%d,%d]\n",proc,i,j);
				}
			}
		}else{
			double *sum = (double*)malloc((end-start)*sizeof(double));
			for (i = start; i < end; i++) {
				double currSum = 0;
				for (k = 0; k < j; k++) {
					currSum = currSum + L[i*n+k] * U_2[k*n+j];	 
					// printf("Mult L[%d,%d](%lf) and U[%d,%d](%lf) for L[%d,%d] with A[%d][%d] by rank %d\n",i,k,L[i*n+k],k,j,U_2[k*n+j],i,j,i,j,rank);
				}
				sum[i-start]=currSum;
			}
			// printf("Thread %d Sending %d size\n",rank,end-start );
			MPI_Send(sum,end-start,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
		}

		//SEND L
		if(rank==0){
			for(int proc = 1; proc< comm_sz; proc++){ 
				double sendL[(n-j)];
				for(i=j;i<n;i++){
					sendL[(i-j)] = L[i*n+j];
				}
				MPI_Send(sendL,(n-j),MPI_DOUBLE,proc,1,MPI_COMM_WORLD);
			}
		}else{
			double sendL[(n-j)];
			Recv(sendL, (n-j), MPI_DOUBLE, 0, 1);
			for(i=j;i<n;i++){
				L[i*n+j]=sendL[(i-j)];
			}
			
		}

		// CALC U
		if(rank==0){
			// printf("Thread %d Receving %d size\n",rank,end-start);
			for (i = start; i < end; i++) {
				double currSum = 0;
				for (k = 0; k < j; k++) {
					currSum = currSum + L[j*n+k] * U_2[k*n+i];	
					// printf("Mult L[%d,%d] and U[%d,%d] for U[%d,%d] with A[%d][%d]\n",j,k,k,i,j,i,j,i);
				}
				U_2[j*n+i] = (A[j*n+i] - currSum) / L[j*n+j];
				// printf("Thread %d wrote in U[%d,%d]\n",rank,i,j);
			}
			for(int proc = 1; proc<comm_sz; proc++){
				int rstart = diff * proc + j;
				int rend = (proc==comm_sz-1)? n : diff * (proc+1) +j;
				double recsum[rend-rstart];
				Recv(recsum, -1, MPI_DOUBLE, proc, 0);

				for (i = rstart; i < rend; i++) {
					U_2[j*n+i] = (A[j*n+i] - recsum[i-rstart]) / L[j*n+j];
					// printf("Thread %d wrote in U[%d,%d]\n",proc,i,j);
				}
			}
		}else{
			double *sum = (double*)malloc((end-start)*sizeof(double));
			for (i = start; i < end; i++) {
				double currSum = 0;
				for (k = 0; k < j; k++) {
					currSum = currSum + L[j*n+k] * U_2[k*n+i];	
					// printf("Mult L[%d,%d](%lf) and U[%d,%d](%lf) for U[%d,%d] with A[%d][%d] by rank %d\n",j,k, L[j*n+k],k,i, U_2[k*n+i],j,i,j,i,rank);

				}
				sum[i-start]=currSum;
			}
			// printf("Thread %d Sending %d size\n",rank,end-start );
			MPI_Send(sum,end-start,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
		}

		//SEND U
		if(rank==0){
			for(int proc = 1; proc< comm_sz; proc++){
				double sendU_2[(n-j)];
				for(i=j;i<n;i++){
					sendU_2[(i-j)] = U_2[j*n+i];
				}
				MPI_Send(sendU_2,(n-j),MPI_DOUBLE,proc,1,MPI_COMM_WORLD);
				// MPI_Send(U_2,n*n,MPI_DOUBLE,proc,1,MPI_COMM_WORLD);
			}
			// write_output(L, U_2, n, comm_sz, 4);
		}else{
			double sendU_2[(n-j)];
			Recv(sendU_2, (n-j), MPI_DOUBLE, 0, 1);
			for(i=j;i<n;i++){
				U_2[j*n+i]=sendU_2[(i-j)];
			}
			// Recv(U_2, n*n, MPI_DOUBLE, 0, 1);
		}
		
		// if(j==1) while(1);
		
	}
	if(rank == 0) {
		// printf("---------------L-----------------\n");
		// print_matrix(L, n);
		// printf("---------------U-----------------\n");
		// print_matrix(U_2, n);
		// printf("---------------LU-----------------\n");
		// matrix_mult(L, U_2, n);
		write_output(L, U_2, n, comm_sz, 4);
	}
	// MPI_Barrier(MPI_COMM_WORLD);
}

int main(int argc, char *argv[]) {
	n = atoi(argv[1]);
	char *input_file = argv[2];
	// comm_sz; // Number of Processes
	int rank; // Process number

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0) {
		init_matrix(n, input_file);
	}
	LU_decomp(n, rank, comm_sz);
	MPI_Finalize();
	return 0;
}