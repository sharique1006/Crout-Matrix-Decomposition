import numpy as np
import sys

A_file = sys.argv[1]
L_file = sys.argv[2]
U_file = sys.argv[3]

def read_from_file(filename):
	# mat = mat.toarray()
	mat=[]
	with open(filename,'r') as f:
		for line in f:
			line1=line.rstrip().split(" ")
			mat.append(list(map(float, line1)))
			# f.write("\n")
	return np.asarray(mat)

A = read_from_file(A_file)
L = read_from_file(L_file)
U = read_from_file(U_file)

print(A)
print(L)
print(U)

print('Is L Lower Triangular? {}'.format(np.allclose(LA, np.tril(LA))))
print('Is U Upper Triangular? {}\n'.format(np.allclose(UA, np.triu(U))))