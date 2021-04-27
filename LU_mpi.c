#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

double *A;
float diff;
int j, n, comm_sz;
char *outfile_L, *outfile_U;


void write_output(char *fname, double* arr, int n ){
	FILE *f = fopen(fname, "w");
	for( int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			fprintf(f, "%0.12f ", arr[n*i + j]);
		}
		fprintf(f, "\n");
	}
	fclose(f);
}
void init_matrix(int n, char *input_file) {
	FILE *f;
	f = fopen(input_file, "r");
	A = (double *)malloc(n*n * sizeof(double ));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (!fscanf(f, "%140lf ", &A[i*n+j])) {
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

void startegy4(int n, int rank, int comm_sz) {

	double *L=(double *)malloc(n*n * sizeof(double));
	double *U=(double *)malloc(n*n * sizeof(double));

	for(int i=0;i<n;i++){
		for(int k=0;k<n;k++){
			L[i*n+k]=0;
			U[i*n+k]=0;
		}
	}

	for(int i = 0; i < n; i++) {
		U[i*n + i] = 1;
	}
	int i,k;
	int start,end;

	// if(rank==0){
	// 	for(int proc = 1; proc< comm_sz; proc++) MPI_Send(U,n*n,MPI_DOUBLE,proc,1,MPI_COMM_WORLD);
	// }else{
	// 	Recv(U, n*n, MPI_DOUBLE, 0, 1);
	// }

	for (j = 0; j < n; j++) {
		diff = (float)(n-j)/(float)(comm_sz);
		start = (diff)*rank+j; 
		end = (rank==comm_sz-1) ? n : (diff)*(rank+1)+j;

		// Calc L
		if(rank==0){
			// printf("=====================================\n");
			// printf("Thread %d Receving %d size\n",rank,end-start);
			for (i = start; i < end; i++) {
				double currSum = 0;
				for (k = 0; k < j; k++) {
					currSum = currSum + L[i*n+k] * U[k*n+j];
					// printf("Mult L[%d,%d](%lf) and U[%d,%d](%lf) for L[%d,%d] with A[%d][%d](%lf)\n",i,k,L[i*n+k],k,j,U[k*n+j],i,j,i,j,A[(i*n+j)]);
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
					currSum = currSum + L[i*n+k] * U[k*n+j];	 
					// printf("Mult L[%d,%d](%lf) and U[%d,%d](%lf) for L[%d,%d] with A[%d][%d] by rank %d\n",i,k,L[i*n+k],k,j,U[k*n+j],i,j,i,j,rank);
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
		if (L[j*n+j] == 0) {	
			printf("Exiting!\n");			
			exit(0);
		}
		// CALC U
		if(rank==0){
			// printf("Thread %d Receving %d size\n",rank,end-start);
			for (i = start; i < end; i++) {
				double currSum = 0;
				for (k = 0; k < j; k++) {
					currSum = currSum + L[j*n+k] * U[k*n+i];	
					// printf("Mult L[%d,%d] and U[%d,%d] for U[%d,%d] with A[%d][%d]\n",j,k,k,i,j,i,j,i);
				}
				U[j*n+i] = (A[j*n+i] - currSum) / L[j*n+j];
				// printf("Thread %d wrote in U[%d,%d]\n",rank,i,j);
			}
			for(int proc = 1; proc<comm_sz; proc++){
				int rstart = diff * proc + j;
				int rend = (proc==comm_sz-1)? n : diff * (proc+1) +j;
				double recsum[rend-rstart];
				Recv(recsum, -1, MPI_DOUBLE, proc, 0);

				for (i = rstart; i < rend; i++) {
					U[j*n+i] = (A[j*n+i] - recsum[i-rstart]) / L[j*n+j];
					// printf("Thread %d wrote in U[%d,%d]\n",proc,i,j);
				}
			}
		}else{
			double *sum = (double*)malloc((end-start)*sizeof(double));
			for (i = start; i < end; i++) {
				double currSum = 0;
				for (k = 0; k < j; k++) {
					currSum = currSum + L[j*n+k] * U[k*n+i];	
					// printf("Mult L[%d,%d](%lf) and U[%d,%d](%lf) for U[%d,%d] with A[%d][%d] by rank %d\n",j,k, L[j*n+k],k,i, U[k*n+i],j,i,j,i,rank);

				}
				sum[i-start]=currSum;
			}
			// printf("Thread %d Sending %d size\n",rank,end-start );
			MPI_Send(sum,end-start,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
		}

		//SEND U
		if(rank==0){
			for(int proc = 1; proc< comm_sz; proc++){
				double sendU[(n-j)];
				for(i=j;i<n;i++){
					sendU[(i-j)] = U[j*n+i];
				}
				MPI_Send(sendU,(n-j),MPI_DOUBLE,proc,1,MPI_COMM_WORLD);
				// MPI_Send(U,n*n,MPI_DOUBLE,proc,1,MPI_COMM_WORLD);
			}
			// write_output(L, U, n, comm_sz, 4);
		}else{
			double sendU[(n-j)];
			Recv(sendU, (n-j), MPI_DOUBLE, 0, 1);
			for(i=j;i<n;i++){
				U[j*n+i]=sendU[(i-j)];
			}
			// Recv(U, n*n, MPI_DOUBLE, 0, 1);
		}
		
		// if(j==1) while(1);
		
	}
	if(rank == 0) {
		write_output(outfile_L, L, n);
		write_output(outfile_U, U, n);
	}
}

int main(int argc, char *argv[]) {
	n = atoi(argv[1]);
	char *input_file = argv[2];
	outfile_L = argv[3];
	outfile_U = argv[4];
	int rank;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0) {
		init_matrix(n, input_file);
	}
	startegy4(n, rank, comm_sz);
	MPI_Finalize();
	return 0;
}