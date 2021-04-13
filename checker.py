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

Adash = np.matmul(L, U)

correct = 1

for i in range(A.shape[0]):
    for j in range(A.shape[0]):
        if abs(Adash[i][j] - A[i][j]) > 1e-3:
            correct = 0
            i = A.shape[0]
            break

if correct == 1:
    print('Correct Decomposition')
else:
    print('Incorrect Decomposition')
    exit(0)

if (abs(np.linalg.det(U) - 1) < 1e-3):
    print('Correct determinant of U')

print('Is L Lower Triangular? {}'.format(np.allclose(L, np.tril(L))))
print('Is U Upper Triangular? {}'.format(np.allclose(U, np.triu(U))))
