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

	def vary_P1_P3(self, min = 0, max = 1, step = 0.05, P2_val = 0.5):
		self.P1_vals = np.arange(min,max,step)
		self.P2_val = P2_val
		self.P3_vals = np.arange(min,max,step)
		self.step = step
		self.I_frac = np.zeros((len(self.P1_vals), len(self.P3_vals)))
		self.I_var = np.zeros((len(self.P1_vals), len(self.P3_vals)))
		for P1 in self.P1_vals:
			for P3 in self.P3_vals:
				print P1, P3

				I_list = np.zeros(self.measurement_len/self.sweeps_per_measurement)

				self.sim = SIRS(self.N,P1,self.P2_val,P3)

				for sweep in range(self.equlibrium_sweeps):
					for i in range(self.sweep):
						self.sim.update()

				for sweep in range(self.measurement_len):
					for i in range(self.sweep):
						self.sim.update()

					if sweep%self.sweeps_per_measurement == 0: #measure every x sweeps
						# I_list.append(self.get_I(self.sim.array))
						# print sweep/self.sweeps_per_measurement
						I_list[sweep/self.sweeps_per_measurement] = self.get_I(self.sim.array)
						# print I_list
				self.I_frac[P1/self.step, P3/self.step] = np.mean(I_list)/self.N
				# print I_list**2
				# var = (np.mean(I_list**2) - np.mean(I_list)**2)/self.N
				# print var
				self.I_var[P1/self.step, P3/self.step] = (np.mean(I_list**2) - np.mean(I_list)**2)/self.N

		np.savetxt(str("Data/Vary_P1_P3_" + str(self.step)  + ".txt"), self.I_frac)
		np.savetxt(str("Data/Vary_P1_P3_Var_" + str(self.step)  + ".txt"), self.I_var)
		return self.I_frac, self.I_var

	def vary_P1(self, min = 0.2, max = 0.5, step = 0.05, P2_val = 0.5, P3_val = 0.2):
		self.P1_vals = np.arange(min,max,step)
		self.P2_val = P2_val
		self.P3_val = P3_val
		self.step = step
		self.I_frac = np.zeros((len(self.P1_vals)))
		# print len(self.I_frac)
		# print self.I_frac
		self.I_var = np.zeros((len(self.P1_vals)))
		step = 0
		for P1 in self.P1_vals:
			# print P1

			I_list = np.zeros(self.measurement_len/self.sweeps_per_measurement)

			self.sim = SIRS(self.N,P1,self.P2_val,self.P3_val)

			for sweep in range(self.equlibrium_sweeps):
				for i in range(self.sweep):
					self.sim.update()

			for sweep in range(self.measurement_len):
				for i in range(self.sweep):
					self.sim.update()

				if sweep%self.sweeps_per_measurement == 0: #measure every x sweeps
					I_list[sweep/self.sweeps_per_measurement] = self.get_I(self.sim.array)

			# print P1/self.step
			print step
			self.I_frac[step] = np.mean(I_list)/self.N
			# var = (np.mean(I_list**2) - np.mean(I_list)**2)/self.N
			# print var
			self.I_var [step] = (np.mean(I_list**2) - np.mean(I_list)**2)/self.N
			step +=1

		np.savetxt(str("Data/Vary_P1_" + str(self.step)  + ".txt"), self.I_frac)
		np.savetxt(str("Data/Vary_P1_Var_" + str(self.step)  + ".txt"), self.I_var)
		return self.I_frac, self.I_var, self.P1_vals

# class Plot():
# 	def __init__(self):
# 		pass
# 	def plot_heatmap(self,data):
# 		plt.imshow(data, interpolation = "nearest")
# 		plt.colorbar()
# 		plt.show()

def main():
	# I_frac, I_var = Experiments(50, 100, 1000, 10).vary_P1_P3(step = 0.2)
	#cut
	I_frac, I_var, P1_vals = Experiments(50, 100, 1000, 10).vary_P1(step = 0.05)
	# Plot().plot_heatmap(data)

main()
