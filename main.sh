#!/bin/bash
bash compile.sh
for matrix_size in {1000..5000..1000}
do
    echo "Matrix Size: $matrix_size"
    for strategy in {1..3..1}
    do
        echo "Strategy: $strategy"
        echo "===================="
        for ((num_threads = 2 ; num_threads <= 16 ; num_threads*=2)); 
        do
            echo "Number of Threads: $num_threads"
            time bash run.sh $matrix_size "input_${matrix_size}.txt" $num_threads $strategy
            python3 format_checker.py "input_${matrix_size}.txt" "output_L_${strategy}_${num_threads}.txt" "output_U_${strategy}_${num_threads}.txt"
            echo ""
        done
        echo "===================="
    done
done


for matrix_size in {1000..5000..1000}
do
    echo "Matrix Size: $matrix_size"
    echo "Strategy: 4"
    echo "===================="
    for ((num_threads = 2 ; num_threads <= 6 ; num_threads+=2)); 
    do
        echo "Number of Threads: $num_threads"
        time bash run.sh $matrix_size "input_${matrix_size}.txt" $num_threads 4
        python3 format_checker.py "input_${matrix_size}.txt" "output_L_4_${num_threads}.txt" "output_U_4_${num_threads}.txt"
        echo ""
    done
    echo "===================="
done
