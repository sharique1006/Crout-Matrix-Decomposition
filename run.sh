make clean
make
make sequential n=50 f="input_50.txt" t=8
python3 checker.py input_50.txt output_L_8_0.txt output_U_8_0.txt
make strategy2 n=50 f="input_50.txt" t=8
python3 checker.py input_50.txt output_L_8_2.txt output_U_8_2.txt
