import numpy as np
import sys

n = int(sys.argv[1])
file = sys.argv[2]

A = np.random.uniform(-4, 4, (n, n))

f = open(file, 'w')
for i in A:
    for j in i:
        print(round(j, 12), '', end='', file=f)
    print(file=f)
f.close()