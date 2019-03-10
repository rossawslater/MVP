import numpy as np
import matplotlib.pyplot as plt


class Plot():
	def __init__(self, file):
		self.data = np.loadtxt(file)

	def plot_heatmap(self):
		plt.imshow(self.data[1:,1:], interpolation = "nearest")
		plt.gca().invert_yaxis()
		plt.colorbar()
		plt.show()

	def plot_graph(self):
		plt.plot(np.arange(0.2,0.5,0.05), self.data)
		plt.show()

def main():
	# Plot("Data/Vary_P1_P3_Var_0.1.txt").plot_heatmap()
	Plot("Data/Vary_P1_0.05.txt").plot_graph()
main()
