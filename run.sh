#!/bin/bash
matrix_size=$1
input_file=$2
num_threads=$3
strategy=$4
outfile_L="output_L_${strategy}_${num_threads}.txt"
outfile_U="output_U_${strategy}_${num_threads}.txt"

if [ ${strategy} == 4 ] ; then
    mpiexec -n $num_threads ./LU_mpi $matrix_size $input_file $outfile_L $outfile_U
else
    ./LU_omp $matrix_size $input_file $num_threads $strategy $outfile_L $outfile_U
fi