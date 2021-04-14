make clean
make
make sequential n=50 f="input_50.txt" t=4
python3 checker.py input_50.txt output_L_4_0.txt output_U_4_0.txt
make strategy1 n=50 f="input_50.txt" t=4
python3 checker.py input_50.txt output_L_4_1.txt output_U_4_1.txt
make strategy2 n=50 f="input_50.txt" t=4
python3 checker.py input_50.txt output_L_4_2.txt output_U_4_2.txt
