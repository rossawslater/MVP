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


def main():
	Plot("Data/Vary_P1_P3_0.1.txt").plot_heatmap()
main()
