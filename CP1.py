import numpy as np
import time
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys
import datetime

class IsingModel():

	def __init__(self, N, T, dynamics, init_state, J = 1):
		self.N = N
		self.T = T
		self.J = J
		self.plotFigure = plt.figure() #required for animation

		if dynamics == "Glauber":
			self.update = self.Glauber
		elif dynamics == "Kawasaki":
			self.update = self.Kawasaki
		else:
			print "Input Error"
			exit()

		if init_state == "Random":
			self.initalize_random()
		elif init_state == "Uniform":
			self.initalize_uniform(dynamics)
		else:
			print "Input Error"
			exit()

		# if type(N) or type(T) != float():
		# 	print "Input Error"
		# 	exit()

	def initalize_uniform(self, dynamics):
		self.array = np.ones([self.N, self.N])
		if dynamics == "Kawasaki":
			for i in range(self.N):
				for j in range(self.N/2):
						self.array[i,j] = -1

	def initalize_random(self):
		self.array = np.zeros([self.N, self.N])
		for x in range(0, self.N):
			for y in range(0, self.N):
				self.array[x,y] = np.random.choice([-1,1])

	def get_x_y(self):
		return np.random.randint(0,self.N), np.random.randint(0,self.N)

	def flip_spin(self, x, y):
		self.array[x,y] *= -1

	def get_dE(self, x, y):

		i = self.array[x, y]

		if x == self.N-1: #boundary conditions so nearest neighbor loops over at far edge
			x_plus_1 = 0
		else:
			x_plus_1 = x + 1

		if y == self.N-1:
			y_plus_1 = 0
		else:
			y_plus_1 = y + 1

		x_minus_1 = x - 1 #if x or y is 0 then -1 indexes to last entry in array so loops over
		y_minus_1 = y - 1

		neighbors = np.array([self.array[x_plus_1, y], self.array[x_minus_1,y], self.array[x, y_plus_1], self.array[x,y_minus_1]])

		nn_sum = np.sum(neighbors)

		E_i = i * nn_sum

		E_f = - i * nn_sum #flip spin of i

		return -E_f + E_i

	def get_probability(self, E):
		r = np.random.rand()
		p = np.exp(-E/self.T)

		if E <= 0:
			return True
		elif r <= p:
			return True
		else:
			return False

	def Glauber(self):
		x,y = self.get_x_y()
		dE = self.get_dE(x,y)
		if self.get_probability(dE) == True:
			self.flip_spin(x,y)
		else:
			pass

	def Kawasaki(self):
		i_x, i_y = self.get_x_y()
		j_x, j_y = self.get_x_y()
		if self.array[i_x,i_y] != self.array[j_x,j_y]:
			dE_1 = self.get_dE(i_x,i_y)
			dE_2 = self.get_dE(j_x,j_y)
			if self.get_probability(dE_1 + dE_2) == True:
				self.flip_spin(i_x,i_y)
				self.flip_spin(j_x,j_y)
		else:
			pass

	def updatePlot(self,i): #move to plotting class?
		for i in range(self.N**2):
			self.update()
		self.plotFigure.clear()# Clear the old plot
		plt.imshow(self.array, interpolation = "nearest", cmap = "binary")# Make the new plot
		plt.axis('off')

	def Visualise(self):# Function that runs the animaion
		ani = animation.FuncAnimation(self.plotFigure, self.updatePlot)
		plt.show()

	# def get_E_no_DC(self, x, y): #move to Experiments class
	# 	i = self.array[x, y]
	#
	# 	if x == self.N-1: #boundary conditions so nearest neighbor loops over at far edge
	# 		x_plus_1 = 0
	# 	else:
	# 		x_plus_1 = x + 1
	#
	# 	if y == self.N-1:
	# 		y_plus_1 = 0
	# 	else:
	# 		y_plus_1 = y + 1
	#
	# 	E = i * np.sum(np.array([self.array[x_plus_1, y], self.array[x, y_plus_1]])) #spin of i * sum of neaest neighbor spins
	#
	# 	return E

class Experiments():

	def __init__(self, N, dynamics, min = 1, max = 3, step = 0.1):
		self.N = N
		self.T_vals = np.arange(min,max,step)
		self.dynamics = dynamics
		self.init_state = "Uniform" #just define in instance of creating object?

	def get_M(self,array):
		return np.sum(array)

	def get_Chi(self,M_Vals,T):
		return (np.mean(M_Vals**2) - np.mean(M_Vals)**2)/(self.N * T)

	def get_total_E(self, array):
		E = 0
		for x in range(self.N):
			for y in range(self.N):
				s = array[x, y]

				if x == self.N-1: #boundary conditions so nearest neighbor loops over at far edge
					x_plus_1 = 0
				else:
					x_plus_1 = x + 1

				if y == self.N-1:
					y_plus_1 = 0
				else:
					y_plus_1 = y + 1

				E += s * np.sum(np.array([array[x_plus_1, y], array[x, y_plus_1]])) #spin of i * sum of neaest neighbor spins

		return E

	def get_C(self,E_Vals,T):
		return (np.mean(E_Vals**2) - np.mean(E_Vals)**2)/(self.N * (T**2))

	def vary_T(self):
		self.Magnetisation = np.zeros((len(self.T_vals), 1000))
		self.Energy = np.zeros((len(self.T_vals), 1000))
		self.Chi = []
		self.C = []

		equlibrium_sweeps = 100
		measurement_sweeps = 10000

		sweep_counter = 0

		step = 0
		for T in self.T_vals:

			self.sim = IsingModel(self.N, T, self.dynamics, "Uniform")

			for sweep in range(equlibrium_sweeps):  #equlibriate
				for i in range(self.N**2):
					self.sim.update()

				sweep_counter += 1
				sys.stdout.write(" Simulation progress: %.2f%%   \r" %(sweep_counter*100/float(len(self.T_vals)*(equlibrium_sweeps+measurement_sweeps))))
				sys.stdout.flush()

			for sweep in range(measurement_sweeps): #run 10000 sweep sample
				for i in range(self.N**2): #sweep is N^2 updates
					self.sim.update()

				if sweep%10 == 0: #measure every 10 sweeps
					self.Magnetisation[step,sweep/10] = self.get_M(self.sim.array)
					self.Energy[step,sweep/10] = self.get_total_E(self.sim.array)

				sweep_counter += 1
				sys.stdout.write(" Simulation progress: %.2f%%   \r" %(sweep_counter*100/float(len(self.T_vals)*(equlibrium_sweeps+measurement_sweeps))))
				sys.stdout.flush()

			self.Chi.append(self.get_Chi(self.Magnetisation[step,:], T))
			self.C.append(self.get_C(self.Energy[step,:], T))

			step += 1 #finished measurements for one value of T ***ONLY used for tracking in terminal output***

		return self.Energy, self.Magnetisation, self.C, self.Chi


	def get_error(self, Data, func): #using bootstrap algorithm ***DECIDE IF ALL self. OR NOT***
		no_steps = len(Data) #give better name
		no_new_measurements = 100 #define what this should be
		no_measurements = len(Data[0])
		new_measurements = np.zeros((no_steps, no_new_measurements))
		new_samples = np.zeros((no_steps, no_measurements))
		stdevs = np.zeros((no_steps))
		# self.T_vals = np.arange(1,3,0.1) #draw this from

		for i in range(no_steps):
			for k in range(no_new_measurements): #define what this should be, 100?
				for j in range(no_measurements):

					new_samples[i] = np.random.choice(Data[i], size = len(Data[0]))
				new_measurements[i,k] = func(new_samples[i], self.T_vals[i])

			stdevs[i] = np.sqrt(np.mean(new_measurements[i]**2) - np.mean(new_measurements[i])**2)

		return stdevs

	def get_mean_vals(self):
		mean_M = []
		for i in range(len(self.Magnetisation)):
			mean_M.append(abs(np.mean(self.Magnetisation[i]))) #absolute value of mean of magnetisation

		mean_E = []
		for i in range(len(self.Energy)):
			mean_E.append(np.mean(self.Energy[i]))

		return mean_E, mean_M

	def get_mean_vals(self, Data):
		means_vals = []
		for i in range(len(Data)):
			mean_vals.append(np.mean(Data[i]))

		return mean_vals

	def output(self):
		pass

def main():
	N = int(sys.argv[1])
	T = float(sys.argv[2])
	dynamics = sys.argv[3]
	init_state = sys.argv[4]

	# IsingModel(N,T,dynamics,init_state).Visualise()
	Energy, Magnetisation, C, Chi =  Experiments(N,dynamics).vary_T()
	#
	# date_time = str(datetime.datetime.now())
	# np.savetxt(str("Data/Energy " + date_time +  ".txt"), Energy)
	# np.savetxt(str("Data/Magnetisation " + date_time +  ".txt"), Magnetisation)
	# np.savetxt(str("Data/Suscepibility " + date_time +  ".txt"), Chi)
	# np.savetxt(str("Data/Heat Capacity " + date_time +  ".txt"), C)
	#
	#
	# np.savetxt("Data/C_Error" + date_time +  ".txt", self.get_error(Energy, x.get_C))
	# np.savetxt("Data/Chi_Error" + date_time +  ".txt", self.get_error(Magnetisation, x.get_Chi))

main()
