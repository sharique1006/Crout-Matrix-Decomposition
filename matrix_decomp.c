#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

double **A, **L, **U;

void displayALU(int n) {
	printf("---------------A-----------------\n");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			printf("%.12lf ", A[i][j]);
		}
		printf("\n");
	}
	printf("---------------L-----------------\n");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			printf("%.12lf ", L[i][j]);
		}
		printf("\n");
	}
	printf("---------------U-----------------\n");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			printf("%.12lf ", U[i][j]);
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
			printf("%.12lf ", sum);
		}
		printf("\n");
	}
}

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
	#pragma omp parallel shared(A,L,U)
	{
	
		for (int i = 0; i < n; i++)
		{
			//for each row....
			//rows are split into seperate threads for processing
			#pragma omp for schedule(static)
			for (int j = 0; j < n; j++)
			{
				//if j is smaller than i, set l[j][i] to 0
				if (j < i)
				{
					L[j][i] = 0;
					continue;
				}
				//otherwise, do some math to get the right value
				L[j][i] = A[j][i];
				for (int k = 0; k < i; k++)
				{
					//deduct from the current l cell the value of these 2 values multiplied
					L[j][i] = L[j][i] - L[j][k] * U[k][i];
				}
			}
			//for each row...
			//rows are split into seperate threads for processing
			#pragma omp for schedule(static)
			for (int j = 0; j < n; j++)
			{
				//if j is smaller than i, set u's current index to 0
				if (j < i)
				{
					U[i][j] = 0;
					continue;
				}
				//if they're equal, set u's current index to 1
				if (j == i)
				{
					U[i][j] = 1;
					continue;
				}
				//otherwise, do some math to get the right value
				U[i][j] = A[i][j] / L[i][i];
				for (int k = 0; k < i; k++)
				{
					U[i][j] = U[i][j] - ((L[i][k] * U[k][j]) / L[i][i]);
				}
			
			}
		}
	}
}

void strategy2(int n, int num_threads) {
	
}

void strategy3(int n, int num_threads) {

}

void strategy4(int n, int num_threads) {

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

void plotGraph(const char gg[], int n) {
    FILE * pipe = popen ("gnuplot -persistent", "w");	
	fprintf(pipe,"set style data lines\n");
	fprintf(pipe, "set title \"Time vs Number of threads, n = %d\"\n", n);
	fprintf(pipe, "set xlabel \"Number of Threads\"\n");
	fprintf(pipe, "set ylabel \"Time Taken\"\n");
    fprintf(pipe,"plot 'time.txt' using 2:xtic(1) title 'Sequential' lt rgb '#406090', '' using 3 title 'Strategy-1' lt rgb '#40FF00', '' using 4 title 'Strategy-2' lt rgb '#440154', '' using 5 title 'Strategy-3' lt rgb '#40FF00', '' using 6 title 'Strategy-4' lt rgb '#40FF00'\n");
	fprintf(pipe,"set term png\n");
	fprintf(pipe,"set output '%s'\n",gg);	
	fprintf(pipe,"replot\n");
}

void plot(int n, int num_threads) {
	double start_seq, end_seq, start_1, end_1, start_2, end_2, start_3, end_3, start_4, end_4;
	int t = 2;
	FILE *temp = fopen("time.txt", "w");
	while (t < 17) {
		start_seq = omp_get_wtime(); sequential(n); end_seq = omp_get_wtime();
		start_1 = omp_get_wtime(); strategy1(n, num_threads); end_1 = omp_get_wtime();
		start_2 = omp_get_wtime(); strategy2(n, num_threads); end_2 = omp_get_wtime();
		start_3 = omp_get_wtime(); strategy3(n, num_threads); end_3 = omp_get_wtime();
		start_4 = omp_get_wtime(); strategy4(n, num_threads); end_4 = omp_get_wtime();
		fprintf(temp, "%d %f %f %f, %f, %f\n", t, (end_seq-start_seq), (end_1-start_1), (end_2-start_2), (end_3-start_3), (end_4-start_4));
		t = t*2;
	}
	plotGraph("plot.png", n);
}

int main(int argc, char *argv[]) {
	int n = atoi(argv[1]);
	char *input_file = argv[2];
	int num_threads = atoi(argv[3]);
	int strategy = atoi(argv[4]);

	init_matrix(n, input_file);
	displayALU(n);

	double start, end;

	if (strategy == 0) { start = omp_get_wtime(); sequential(n); end = omp_get_wtime(); }
	else if (strategy == 1) { start = omp_get_wtime(); strategy1(n, num_threads); end = omp_get_wtime(); }
	else if (strategy == 2) { start = omp_get_wtime(); strategy2(n, num_threads); end = omp_get_wtime(); }
	else if (strategy == 3) { start = omp_get_wtime(); strategy3(n, num_threads); end = omp_get_wtime(); }
	else { start = omp_get_wtime(); strategy4(n, num_threads); end = omp_get_wtime(); }

	displayLU(n);

	double time_taken = (end - start);
	printf("\nTime taken in strategy %d is %f seconds\n", strategy, time_taken);
	
	write_to_file(n, num_threads, strategy);
	//plot(n, num_threads);
	return 0;
}