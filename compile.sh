#!/bin/bash
gcc -o0 -fopenmp LU_omp.c -o LU_omp
mpicc -g -Wall -o arrLU_mpi arrLU_mpi.c