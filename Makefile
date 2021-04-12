all: compile

compile:
	@ gcc -o0 -fopenmp matrix_decomp.c -o matrix_decomp

sequential:
	@ ./matrix_decomp 4 input.txt 8 0

strategy1:
	@ ./matrix_decomp 4 input.txt 8 1 

strategy2:
	@ ./matrix_decomp 4 input.txt 8 2

strategy3:
	@ ./matrix_decomp 4 input.txt 8 3 

strategy4:
	@ ./matrix_decomp 4 input.txt 8 4 

clean:
	@ rm -rf output*
	@ rm -rf input_*
	@ rm -rf matrix_decomp