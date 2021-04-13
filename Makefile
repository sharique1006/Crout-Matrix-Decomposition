all: compile

compile:
	@ gcc -o0 -fopenmp matrix_decomp.c -o matrix_decomp

sequential:
<<<<<<< HEAD
	@ ./matrix_decomp 500 input_500.txt 8 0

strategy1:
	@ ./matrix_decomp 500 input_500.txt 8 1 

strategy2:
	@ ./matrix_decomp 500 input_500.txt 8 2

strategy3:
	@ ./matrix_decomp 500 input_500.txt 8 3 

strategy4:
	@ ./matrix_decomp 500 input_500.txt 8 4 
=======
	@ ./matrix_decomp $(n) $(f) $(t) 0

strategy1:
	@ ./matrix_decomp $(n) $(f) $(t) 1 

strategy2:
	@ ./matrix_decomp $(n) $(f) $(t) 2

strategy3:
	@ ./matrix_decomp $(n) $(f) $(t) 3 

strategy4:
	@ ./matrix_decomp $(n) $(f) $(t) 4 
>>>>>>> 64c01104fc8c4365324074def3649c79af672426

clean:
	@ rm -rf output*
	@ rm -rf input_*
	@ rm -rf time.txt
	@ rm -rf plot.png
	@ rm -rf matrix_decomp