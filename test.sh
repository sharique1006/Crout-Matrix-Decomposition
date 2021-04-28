#!/bin/bash
make clean
N=1000
num_threads=7
input_file="input_${N}.txt"
python3 gen_matrix.py $N $input_file

echo Data Generated

strategy=1
outfile_L="output_L_${strategy}_${num_threads}.txt"
outfile_U="output_U_${strategy}_${num_threads}.txt"
bash compile.sh
bash run.sh $N $input_file $num_threads $strategy
python3 format_checker.py $input_file $outfile_L $outfile_U

strategy=2
outfile_L="output_L_${strategy}_${num_threads}.txt"
outfile_U="output_U_${strategy}_${num_threads}.txt"
bash compile.sh
bash run.sh $N $input_file $num_threads $strategy
python3 format_checker.py $input_file $outfile_L $outfile_U

strategy=3
outfile_L="output_L_${strategy}_${num_threads}.txt"
outfile_U="output_U_${strategy}_${num_threads}.txt"
bash compile.sh
bash run.sh $N $input_file $num_threads $strategy
python3 format_checker.py $input_file $outfile_L $outfile_U

strategy=4
outfile_L="output_L_${strategy}_${num_threads}.txt"
outfile_U="output_U_${strategy}_${num_threads}.txt"
bash compile.sh
time bash run.sh $N $input_file $num_threads $strategy
python3 format_checker.py $input_file $outfile_L $outfile_U