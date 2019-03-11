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
		# self.sweep = self.N **2

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
				print (P1, P3)

				I_list = np.zeros(self.measurement_len/self.sweeps_per_measurement)

				self.sim = SIRS(self.N,P1,self.P2_val,P3)

				for sweep in range(self.equlibrium_sweeps):
					# for i in range(self.sweep):
					self.sim.update()

				for sweep in range(self.measurement_len):
					# for i in range(self.sweep):
					self.sim.update()

					if sweep%self.sweeps_per_measurement == 0: #measure every x sweeps

						I_list[sweep/self.sweeps_per_measurement] = self.get_I(self.sim.array)

				self.I_frac[P1/self.step, P3/self.step] = float(np.mean(I_list))/self.N**2
				print(self.I_frac[P1/self.step, P3/self.step])
				var = float((np.mean(I_list**2) - np.mean(I_list)**2))/self.N**2
				print (var)
				self.I_var[P1/self.step, P3/self.step] = var# (np.mean(I_list**2) - np.mean(I_list)**2)/self.N**2

		print (self.I_var)

		np.savetxt(str("Data/Vary_P1_P3_" + str(self.step)  + ".txt"), self.I_frac)
		np.savetxt(str("Data/Vary_P1_P3_Var_" + str(self.step)  + ".txt"), self.I_var)
		return self.I_frac, self.I_var, self.P1_vals, self.P3_vals

	def vary_P1(self, min = 0.2, max = 0.5, step = 0.05, P2_val = 0.5, P3_val = 0.2):
		self.P1_vals = np.arange(min,max,step)
		self.P2_val = P2_val
		self.P3_val = P3_val
		self.step = step
		self.I_frac = np.zeros((len(self.P1_vals)))
		self.I_var = np.zeros((len(self.P1_vals)))

		step = 0

		for P1 in self.P1_vals:
			print (P1)

			I_list = np.zeros(self.measurement_len/self.sweeps_per_measurement)

			self.sim = SIRS(self.N,P1,self.P2_val,self.P3_val)

			for sweep in range(self.equlibrium_sweeps):
				# for i in range(self.sweep):
				self.sim.update()

			for sweep in range(self.measurement_len):
				# for i in range(self.sweep):
				self.sim.update()

				if sweep%self.sweeps_per_measurement == 0: #measure every x sweeps
					I_list[sweep/self.sweeps_per_measurement] = self.get_I(self.sim.array)

			self.I_frac[step] = float(np.mean(I_list))/self.N**2
			self.I_var[step] = float(np.mean(I_list**2) - np.mean(I_list)**2)/self.N**2

			step +=1

		np.savetxt(str("Data/Vary_P1_" + str(self.step)  + ".txt"), self.I_frac)
		np.savetxt(str("Data/Vary_P1_Var_" + str(self.step)  + ".txt"), self.I_var)

		return self.I_frac, self.I_var, self.P1_vals

	def vary_immunity(self, min = 0, max = 1, step = 0.05):
		immunity_fracs = np.arange(min,max,step)
		self.I_frac = np.zeros((len(immunity_fracs)))
		self.Error = np.zeros((len(immunity_fracs)))
		self.step = step
		step = 0
		for frac in immunity_fracs:
			print (frac)
			I_list = np.zeros(self.measurement_len/self.sweeps_per_measurement)

			self.sim = SIRS(self.N,0.5,0.5,0.5,frac)

			for sweep in range(self.equlibrium_sweeps):
				self.sim.update()

			for sweep in range(self.measurement_len):
				self.sim.update()

				if sweep%self.sweeps_per_measurement == 0: #measure every x sweeps
					print (float(self.get_I(self.sim.array))/self.N**2)
					I_list[sweep/self.sweeps_per_measurement] = self.get_I(self.sim.array)

			self.I_frac[step] = float(np.mean(I_list))/self.N**2
			self.Error[step] = np.std(I_list)/np.sqrt(len(I_list))

			step += 1

		np.savetxt(str("Data/Vary_Immune_" + str(self.step)  + ".txt"), self.I_frac)
		np.savetxt(str("Data/Vary_Immune_" + str(self.step)  + "_errors.txt"), self.Error)

		return self.I_frac, immunity_fracs, self.Error


# class Plot():
# 	def __init__(self):
# 		pass
# 	def plot_heatmap(self,data):
# 		plt.imshow(data, interpolation = "nearest")
# 		plt.colorbar()
# 		plt.show()

def main():
	I_frac, I_var, P1_vals, P3_vals = Experiments().vary_P1_P3(step = 0.05)
	#cut
	# I_frac, I_var, P1_vals = Experiments().vary_P1(step = 0.05)
	# Plot().plot_heatmap(data)
	# I_frac, immunity_fracs, Error = Experiments().vary_immunity()

main()
