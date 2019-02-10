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

	def initalize_uniform(self, dynamics):
		# if dynamics == "Glauber":
		self.array = np.ones([self.N, self.N])
		if dynamics == "Kawasaki":
			self.array = np.ones([self.N, self.N])
			for i in range(self.N):
				for j in range(self.N):
					if j <= self.N/2:
						self.array[i,j] = 1#np.ones([self.N, self.N]) #change to half and half for kawasaki?
					else:
						self.array[i,j] = -1 #np.ones([self.N, self.N]) #change to half and half for kawasaki?

	def initalize_random(self):
		self.array = np.zeros([self.N, self.N])
		for x in range(0, self.N):
			for y in range(0, self.N):
				self.array[x,y] = np.random.choice([-1,1])

	def get_x_y(self):
		return np.random.randint(0,self.N), np.random.randint(0,self.N)

	def flip_spin(self, x, y):#, x, y):
		self.array[x,y] *= -1

	def get_flip_E(self, x, y):

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

		return E_i, E_f

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
		E_i, E_f = self.get_flip_E(x,y)
		if self.get_probability(-E_f+E_i) == True:
			self.flip_spin(x,y)
		else:
			pass

	def Kawasaki(self):

		i_x, i_y = self.get_x_y()
		j_x, j_y = self.get_x_y()

		if self.array[i_x,i_y] != self.array[j_x,j_y]:
			E_1_init, E_1_final = self.get_flip_E(i_x,i_y)
			E_2_init, E_2_final = self.get_flip_E(j_x,j_y)

			if self.get_probability(-E_1_final + E_1_init + -E_2_final + E_2_init) == True:
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

	def get_E_no_DC(self, x, y): #move to Experiments class

		i = self.array[x, y]

		if x == self.N-1: #boundary conditions so nearest neighbor loops over at far edge
			x_plus_1 = 0
		else:
			x_plus_1 = x + 1

		if y == self.N-1:
			y_plus_1 = 0
		else:
			y_plus_1 = y + 1

		E = i * np.sum(np.array([self.array[x_plus_1, y], self.array[x, y_plus_1]])) #spin of i * sum of neaest neighbor spins

		return E

class Experiments():

	def __init__(self, N, dynamics, min = 1, max = 3):
		self.N = N
		self.T_vals = np.arange(min,max,0.1)
		self.Magnetisation = []

		self.Energy = []
		self.dynamics = dynamics
		self.init_state = "Uniform" #just define in instance of creating object?
		self.Chi = []
		self.C = []

	def get_M(self,array):
		return np.sum(array)

	def get_Chi(self,M_Vals,T):
		return (np.mean(M_Vals**2) - np.mean(M_Vals)**2)/(self.N * T)

	def get_total_E(self):
		E = 0
		for i in range(self.N):
			for j in range(self.N):
				E += self.sim.get_E_no_DC(i,j)

		return E #need function that only calculates E for two neaest neighbours

	def get_C(self,E_Vals,T):
		return (np.mean(E_Vals**2) - np.mean(E_Vals)**2)/(self.N * (T**2))

	def vary_T(self):
		step = 0
		self.Magnetisation = np.zeros((len(self.T_vals), 1000))
		self.Energy = np.zeros((len(self.T_vals), 1000))

		for T in self.T_vals:

			self.sim = IsingModel(self.N, T, self.dynamics, "Uniform")

			for sweep in range(100):  #equlibriate
				for i in range(self.N**2):
					self.sim.update()

			for sweep in range(10000): #run 10000 sweep sample
				for i in range(self.N**2): #sweep is N^2 updates
					self.sim.update()
				if sweep%10 == 0: #measure every 10 sweeps
					self.Magnetisation[step,sweep/10] = self.get_M(self.sim.array)
					self.Energy[step,sweep/10] = self.get_total_E() #need function that gets only one horizontal and one vertical nn Energy

			self.Chi.append(self.get_Chi(self.Magnetisation[step,:], T))
			self.C.append(self.get_C(self.Energy[step,:], T))
			print step
			step += 1 #finished measurements for one value of T

			# sys.stdout.write(" Simulation progress: %.1f%%   \r" %(t*100/float(self.sim_length)))
			# sys.stdout.flush()
		return self.Energy, self.Magnetisation, self.C, self.Chi


	def get_error(self, Data, func): #using bootstrap algorithm
		no_steps = len(Data) #give better name
		no_new_measurements = 100 #define what this should be
		no_measurements = len(Data[0])
		new_measurements = np.zeros((no_steps, no_new_measurements))
		new_samples = np.zeros((no_steps, no_measurements))
		stdevs = np.zeros((no_steps))
		self.T_vals = np.arange(1,3,0.1) #draw this from

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

def main():
	N = int(sys.argv[1])
	T = float(sys.argv[2])
	dynamics = sys.argv[3]
	init_state = sys.argv[4]

	IsingModel(N,T,dynamics,init_state).Visualise()
	# x = Experiments(N,dynamics)
	#
	#
	# Energy, Magnetisation, C, Chi = x.vary_T()
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
