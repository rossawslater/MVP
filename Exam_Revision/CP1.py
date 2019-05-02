import numpy as np
import time
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys


class IsingModel():

	def __init__(self, N, T, dynamics, init_state, J = 1, Kb = 1):
		self.N = int(N)
		self.T = float(T)
		self.J = J #J and Kb are set to 1 unless otherwise stated
		self.Kb = Kb
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
		self.array = np.ones([self.N, self.N])
		if dynamics == "Kawasaki": #need to split array to half 1, half -1 as if all same spin state then no swaps occur
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

	def check_nn(self,i_x,i_y,j_x,j_y): #funcion to check if two array elements are nearest neighbors
		if abs(i_x-j_x) == 1 or abs(i_y-j_y) == 1 or (i_x-j_x)%self.N-1 == 0 or (i_y-j_y)%self.N-1 == 0:
			return True

		else:
			return False

	def get_dE(self, x, y): #calculate the change in energy when the spin is flipped

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

		E_i = i * nn_sum * self.J

		E_f = - i * nn_sum * self.J #flip spin of i

		return -E_f + E_i

	def get_probability(self, E): #using metropolis algorithm returns True or False for flip
		r = np.random.rand()
		p = np.exp(-E/(self.Kb * self.T))

		if E <= 0:
			return True
		elif r <= p:
			return True
		else:
			return False

	def Glauber(self): #glauber dynamics for ising model
		x,y = self.get_x_y()
		dE = self.get_dE(x,y)

		if self.get_probability(dE) == True:
			self.flip_spin(x,y)
		else:
			pass

	def Kawasaki(self): #kawasaki dynamics for ising model
		i_x, i_y = self.get_x_y() #choose two random points
		j_x, j_y = self.get_x_y()
		if self.array[i_x,i_y] != self.array[j_x,j_y]: #if they are not the same
			dE_1 = self.get_dE(i_x,i_y) #find change in energy from flipping spin at each point
			dE_2 = self.get_dE(j_x,j_y)

			if self.check_nn(i_x,i_y,j_x,j_y) == True: # if they are nearest neighbors
				if self.get_probability(dE_1 + dE_2 + 4) == True: #add +4 for double counting
						self.flip_spin(i_x,i_y)#if metropolis algorithm returns true, flip spins
						self.flip_spin(j_x,j_y)
			else:
				if self.get_probability(dE_1 + dE_2) == True:
						self.flip_spin(i_x,i_y)
						self.flip_spin(j_x,j_y)

	def updatePlot(self,i): #function updates plot every 10 sweeps, as like measurements
		for sweep in range(10):
			for i in range(self.N**2):
				self.update()
		self.plotFigure.clear()# Clear the old plot
		plt.imshow(self.array, interpolation = "nearest", cmap = "binary")# Make the new plot
		plt.axis('off')

	def Visualise(self):# Function that runs the animaion
		ani = animation.FuncAnimation(self.plotFigure, self.updatePlot)
		plt.show()

class Experiments():

	def __init__(self, N, dynamics, min = 1, max = 3, step = 0.1):
		self.N = N
		self.T_vals = np.arange(min,max,step)
		self.dynamics = dynamics

	def get_M(self,array): #sum of all spins in array
		return abs(np.sum(array))

	def get_Chi(self,M_Vals,T): #returns suscepibility for given magnetisation measurements and temp
		return (np.mean(M_Vals**2) - np.mean(M_Vals)**2)/(self.N * T)

	def get_total_E(self, array): #find total energy of the array, not double counting interactions
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

	def get_C(self,E_Vals,T): #returns heat capacity for a given energy values and temp
		return (np.mean(E_Vals**2) - np.mean(E_Vals)**2)/(self.N * (T**2))

	def vary_T(self): #experiment to measure different variables as function of T
		Magnetisation = np.zeros((len(self.T_vals), 1000))
		Energy = np.zeros((len(self.T_vals), 1000))
		Chi = np.zeros((len(self.T_vals)))
		C = np.zeros((len(self.T_vals)))

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
					Magnetisation[step,sweep/10] = self.get_M(self.sim.array)
					Energy[step,sweep/10] = self.get_total_E(self.sim.array)

				sweep_counter += 1
				sys.stdout.write(" Simulation progress: %.2f%%   \r" %(sweep_counter*100/float(len(self.T_vals)*(equlibrium_sweeps+measurement_sweeps))))
				sys.stdout.flush()


			Chi[step] = self.get_Chi(Magnetisation[step,:], T)
			C[step] = self.get_C(Energy[step,:], T)

			step += 1 #finished measurements for one value of T ***ONLY used for tracking in terminal output***

		return Energy, Magnetisation, C, Chi


	def get_error(self, Data, func):
		#calculate error on function using bootstrap algorithm
		#must pass in function for use in value calculation
		no_steps = len(Data)
		no_measurements = len(Data[0]) #number of measurements previously made at each T
		no_new_measurements = 100 #take sample of 100 randomly picked samples

		new_samples = np.zeros((no_steps, no_measurements))
		new_measurements = np.zeros((no_steps, no_new_measurements))
		stdevs = np.zeros((no_steps))

		for i in range(no_steps):
			for k in range(no_new_measurements):
				new_samples[i] = np.random.choice(Data[i], size = len(Data[0]))
			new_measurements[i,k] = func(new_samples[i], self.T_vals[i])

			stdevs[i] = np.sqrt(np.mean(new_measurements[i]**2) - np.mean(new_measurements[i])**2)

		return stdevs

	def get_mean_vals(self, Data): #function to calculate mean values of a data set
		mean_vals = []
		for i in range(len(Data)):
			mean_vals.append(np.mean(Data[i]))

		return mean_vals


def main():
	N = int(sys.argv[1]) #arguments passed in in terminal
	T = float(sys.argv[2])
	dynamics = sys.argv[3]
	init_state = sys.argv[4]

	if len(sys.argv) == 6: #option to visualise the model
		visualise = sys.argv[5]

		if visualise == "Visualise":
			IsingModel(N,T,dynamics,init_state).Visualise()
		else:
			print "Input Error"
			exit()
	else:

		exp = Experiments(N,dynamics)
		Energy, Magnetisation, C, Chi =  exp.vary_T()

		np.savetxt(str("Data/Energy_" + str(dynamics)  + ".txt"), Energy)
		np.savetxt(str("Data/Magnetisation_" + str(dynamics) +  ".txt"), Magnetisation)
		np.savetxt(str("Data/Suscepibility_" + str(dynamics) +  ".txt"), Chi)
		np.savetxt(str("Data/Heat_Capacity_" + str(dynamics) + ".txt"), C)

		np.savetxt(str("Data/Mean_Energy_" + str(dynamics)  + ".txt"), exp.get_mean_vals(Energy))
		np.savetxt(str("Data/Mean_Magnetisation_" + str(dynamics)  + ".txt"), exp.get_mean_vals(Magnetisation))

		np.savetxt("Data/C_Error_"+ str(dynamics) + "txt", exp.get_error(Energy, exp.get_C))
		np.savetxt("Data/Chi_Error_" + str(dynamics) + ".txt", exp.get_error(Magnetisation, exp.get_Chi))

main()
