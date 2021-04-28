import matplotlib.pyplot as plt

def plot_seq(T, N):
	plt.plot(N, T_seq)
	plt.title('Time vs N for Sequential LU Decomposition')
	plt.xlabel('Matrix Size')
	plt.ylabel('Run Time in seconds')
	plt.savefig('plot_0.png')
	plt.show()
	plt.close()

T_seq = [2.467, 23.496, 90.096, 230.549, 505.011]
N = [1000, 2000, 3000, 4000, 5000]
plot_seq(T_seq, N)

####### Parallel #######

T1 = [2, 4, 8, 16]
T2 = [2, 4, 6]



N1 = 1000
Strategy1_N1 = [1.657, 1.341, 1.345, 1.531]
Strategy2_N1 = [1.641, 1.633, 2.152, 2.069]
Strategy3_N1 = [1.620, 1.618, 2.045, 2.038]
Strategy4_N1 = [2.343, 1.812, 3.268]

N2 = 2000
Strategy1_N2 = [12.943, 8.933, 10.828, 11.849]
Strategy2_N2 = [13.819, 13.535, 21.548, 18.054]
Strategy3_N2 = [13.188, 15.019, 21.018, 17.898]
Strategy4_N2 = [18.985, 14.408, 12.790]

N3 = 3000
Strategy1_N3 = [46.153, 28.888, 37.945, 30.085]
Strategy2_N3 = [59.510, 89.992, 113.025, 90.496]
Strategy3_N3 = [61.524, 64.513, 68.372, 63.935]
Strategy4_N3 = [54.929, 40.202, 59.603]

N4 = 4000
Strategy1_N4 = [133.853, 72.189, 77.144, 79.706]
Strategy2_N4 = [155.631, 180.614, 214.678, 179.242]
Strategy3_N4 = [132.167, 124.886, 181.654, 167.351]
Strategy4_N4 = [107.574, 85.620, 113.970]

N5 = 5000
Strategy1_N5 = [266.734, 197.782, 199.121, 196.847]
Strategy2_N5 = [263.660, 280.673, 411.029, 414.521]
Strategy3_N5 = [257.201, 260.340, 305.780, 325.557]
Strategy4_N5 = [224.259, 155.35, 230.167]

# Strategy 1, 2, 3
def plot(s, t, s_N1, s_N2, s_N3, s_N4, s_N5):
	plt.plot(t, s_N1, label='N=1000')
	plt.plot(t, s_N2, label='N=2000')
	plt.plot(t, s_N3, label='N=3000')
	plt.plot(t, s_N4, label='N=4000')
	plt.plot(t, s_N5, label='N=5000')
	title = "Run Time vs Number of Threads for Strategy-{}".format(s)
	plt.title(title)
	plt.xlabel('Number of Threads')
	plt.ylabel('Run Time in seconds')
	plt.legend()
	file = "plot_{}.png".format(s)
	plt.savefig(file)
	plt.show()
	plt.close()

# Strategy 4
plot(1, T1, Strategy1_N1, Strategy1_N2, Strategy1_N3, Strategy1_N4, Strategy1_N5)
plot(2, T1, Strategy2_N1, Strategy2_N2, Strategy2_N3, Strategy2_N4, Strategy2_N5)
plot(3, T1, Strategy3_N1, Strategy3_N2, Strategy3_N3, Strategy3_N4, Strategy3_N5)
plot(4, T2, Strategy4_N1, Strategy4_N2, Strategy4_N3, Strategy4_N4, Strategy4_N5)
