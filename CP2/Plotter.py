import numpy as np
import matplotlib.pyplot as plt


class Plot():
	def __init__(self):
		pass

	def read_data(self, file):
		data = np.loadtxt(file)
		return data

	def plot_heatmap(self, data):
		plt.imshow(data, interpolation = "nearest", extent = [0,1,1,0])
		plt.xlabel("P1")
		plt.ylabel("P3")
		plt.title("P1 vs P3 Heatmap")
		plt.gca().invert_yaxis()
		plt.colorbar()
		plt.show()

	def plot_graph(self, data, errors = None):
		plt.errorbar(np.arange(0,1,0.05), data, yerr = errors, xerr = None)
		plt.show()

def main():
	# Plot("Data/Vary_P1_P3_Var_0.05.txt").plot_heatmap()
	p = Plot()
	# data = p.read_data("Data/Vary_Immune_0.05.txt")
	# errors = p.read_data("Data/Vary_Immune_0.05_errors.txt")/2500
	# p.plot_graph(data, errors)


	data = p.read_data("Data/Vary_P1_P3_Var_0.05.txt")
	p.plot_heatmap(data)
main()
