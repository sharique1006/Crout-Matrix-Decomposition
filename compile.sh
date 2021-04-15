#!/bin/bash
gcc -o0 -fopenmp LU_omp.c -o LU_omp
mpicc -g -Wall -o LU_mpi LU_mpi.c