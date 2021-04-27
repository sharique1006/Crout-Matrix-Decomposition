#!/bin/bash
matrix_size=$1
input_file=$2
num_threads=$3
strategy=$4

if [ ${strategy} == 4 ] ; then
    mpiexec -n $num_threads ./arrLU_mpi $matrix_size $input_file
else
    ./LU_omp $matrix_size $input_file $num_threads $strategy
fi