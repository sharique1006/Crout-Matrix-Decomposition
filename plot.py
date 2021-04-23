import matplotlib.pyplot as plt

t = [2, 4, 8, 16]

N1 = 1000
Strategy1_N1 = [1.637, 1.388, 1.416, 1.621]
Strategy2_N1 = [1.637, 1.647, 2.195, 2.242]
Strategy3_N1 = []
Strategy4_N1 = []

N2 = 2000
Strategy1_N2 = [14.648, 9.501, 10.983, 12.098]
Strategy2_N2 = [15.753, 13.930, 20.960, 19.599]
Strategy3_N2 = []
Strategy4_N2 = []

N3 = 3000
Strategy1_N3 = [52.342, 32.802, 41.153, 34.280]
Strategy2_N3 = [50.680, 49.703, 83.528, 76.223]
Strategy3_N3 = []
Strategy4_N3 = []

N4 = 4000
Strategy1_N4 = [133.853, 72.189, 77.144, 79.706]
Strategy2_N4 = [132.167, 124.886, 181.654, 176.351]
Strategy3_N4 = []
Strategy4_N4 = []

N5 = 5000
Strategy1_N5 = []
Strategy2_N5 = []
Strategy3_N5 = []
Strategy4_N5 = []

def plot(N, s1, s2, s3, s4):
	plt.plot(t, s1, label='Strategy-1')
	plt.plot(t, s2, label='Strategy-2')
	plt.plot(t, s3, label='Strategy-3')
	plt.plot(t, s4, label='Strategy-4')
	title = "Run Time vs Number of Threads for N = {}".format(N)
	plt.title(title)
	plt.xlabel('Number of Threads')
	plt.ylabel('Run Time in seconds')
	plt.legend()
	file = "plot_{}.png".format(N)
	plt.savefig(file)
	plt.show()
	plt.close()

plot(N1, Strategy1_N1, Strategy2_N1, Strategy3_N1, Strategy4_N1)
plot(N2, Strategy1_N2, Strategy2_N2, Strategy3_N2, Strategy4_N2)
plot(N3, Strategy1_N3, Strategy2_N3, Strategy3_N3, Strategy4_N3)
plot(N4, Strategy1_N4, Strategy2_N4, Strategy3_N4, Strategy4_N4)
plot(N5, Strategy1_N5, Strategy2_N5, Strategy3_N5, Strategy4_N5)
