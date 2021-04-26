all: omp mpi arr_mpi

omp:
	@ gcc -o0 -fopenmp LU_omp.c -o LU_omp

mpi:
	@ mpicc -g -Wall -o LU_mpi LU_mpi.c

arr_mpi:
	@ mpicc -g -Wall -o arrLU_mpi arrLU_mpi.c
	
sequential:
	@ ./LU_omp $(n) $(f) $(t) 0

strategy1:
	@ ./LU_omp $(n) $(f) $(t) 1 

strategy2:
	@ ./LU_omp $(n) $(f) $(t) 2

strategy3:
	@ ./LU_omp $(n) $(f) $(t) 3 

strategy4:
	@ mpiexec -n $(t) ./LU_mpi $(n) $(f)

arrstrategy4:
	@ mpiexec -n $(t) ./arrLU_mpi $(n) $(f)

clean:
	@ rm -rf LU_omp
	@ rm -rf LU_mpi
	@ rm -rf arrLU_mpi
	@ rm -rf output*