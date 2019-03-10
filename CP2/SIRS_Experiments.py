from CP2 import SIRS
from CP2 import Array
import numpy as np
import matplotlib.pyplot as plt

class Experiments():

	def __init__(self, N = 50, equlibrium_sweeps = 100, measurement_len = 1000, sweeps_per_measurement = 10):
		self.N = N
		self.equlibrium_sweeps = equlibrium_sweeps
		self.measurement_len = measurement_len
		self.sweeps_per_measurement = sweeps_per_measurement
		self.sweep = self.N **2

	def get_I(self,array): #sum infected sites in array
		return np.sum(array == 1)

	def vary_P1_P3(self, min = 0, max = 1, step = 0.05):
		self.P1_vals = np.arange(min,max,step)
		self.P3_vals = np.arange(min,max,step)
		self.P2_val = 0.5
		self.step = step
		self.I_frac = np.zeros((len(self.P1_vals)+1, len(self.P3_vals)+1))
		self.I_var = np.zeros((len(self.P1_vals)+1, len(self.P3_vals)+1))
		for P1 in self.P1_vals:
			for P3 in self.P3_vals:
				print P1, P3

				I_list = []#np.zeros(self.measurement_len/self.measurement_sweeps)

				self.sim = SIRS(self.N,P1,0.5,P3)

				for sweep in range(self.equlibrium_sweeps):
					for i in range(self.sweep):
						self.sim.update()

				for sweep in range(self.measurement_len):
					for i in range(self.sweep):
						self.sim.update()

					if sweep%self.sweeps_per_measurement == 0: #measure every x sweeps
						I_list.append(self.get_I(self.sim.array))

				self.I_frac[P1/self.step +1, P3/self.step + 1] = np.mean(I_list)/self.N
				self.I_frac[P1/self.step, 0] = P1
				self.I_frac[0,P3/self.step] = P3
				# self.I_var [P1/self.step, P3/self.step] = (np.mean(I_list*2) - np.mean(I_list)**2)/self.N

		np.savetxt(str("Data/Vary_P1_P3_" + str(self.step)  + ".txt"), self.I_frac)
		return self.I_frac

	def vary_P1(self, min = 0.2, max = 0.5, step = 0.05):
		self.P1_vals = np.arange(min,max,step)
		self.P3_vals = np.arange(min,max,step)
		self.P2_val = 0.5
		self.step = step
		self.I_frac = np.zeros((len(self.P1_vals), len(self.P3_vals)))
		self.I_var = np.zeros((len(self.P1_vals), len(self.P3_vals)))
		for P1 in self.P1_vals:
			for P3 in self.P3_vals:
				print P1, P3

				I_list = []#np.zeros(self.measurement_len/self.measurement_sweeps)

				self.sim = SIRS(self.N,P1,0.5,P3)

				for sweep in range(self.equlibrium_sweeps):
					for i in range(self.sweep):
						self.sim.update()

				for sweep in range(self.measurement_len):
					for i in range(self.sweep):
						self.sim.update()

					if sweep%self.sweeps_per_measurement == 0: #measure every x sweeps
						I_list.append(self.get_I(self.sim.array))

				self.I_frac[P1/self.step, P3/self.step] = np.mean(I_list)/self.N
				self.I_var [P1/self.step, P3/self.step] = (np.mean(I_list*2) - np.mean(I_list)**2)/self.N

		np.savetxt(str("Data/Vary_P1_P3" + str(self.step)  + ".txt"), self.I_frac)
		return self.I_frac

class Plot():
	def __init__(self):
		pass
	def plot_heatmap(self,data):
		plt.imshow(data, interpolation = "nearest")
		plt.colorbar()
		plt.show()

def main():
	data = Experiments(10, 10, 10, 1).vary_P1_P3(step = 0.1)
	# Plot().plot_heatmap(data)

main()