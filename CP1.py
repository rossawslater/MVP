import numpy as np
import time
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys

class IsingModel():

	def __init__(self, N, T, dynamics, init_state, J = 1):
		self.N = N
		self.T = T
		self.J = J
		self.plotFigure = plt.figure()

		if init_state == "Random":
			self.initalize_random()
		elif init_state == "Uniform":
			self.initalize_uniform()
		else:
			print "Input Error"
			exit()

		if dynamics == "Glauber":
			self.update = self.Glauber
		elif dynamics == "Kawasaki":
			self.update = self.Kawasaki
		else:
			print "Input Error"
			exit()


	def initalize_uniform(self):
		self.array = np.ones([self.N, self.N])
		self.test_array = np.copy(self.array)

	def initalize_random(self):
		self.array = np.zeros([self.N, self.N])
		for x in range(0, self.N):
			for y in range(0, self.N):
				self.array[x,y] = np.random.choice([-1,1])
		self.test_array = np.copy(self.array)

	def get_x_y(self):
		# i = np.random.randint(0,self.N,(2))
		# self.x = i[0]; self.y = i[1]
		# return self.x, self.y
		return np.random.randint(0,self.N), np.random.randint(0,self.N)

	def get_i_and_j(self):  ################## not self?  # tidy so no using self  shouldn't need, just call get i twoice
		r = np.random.randint(0,self.N,(4))
		self.i_x = r[0]; self.i_y = r[1]
		self.j_x = r[2]; self.j_y = r[3]
		return r[0], r[1], r[2], r[3]

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

		i *= -1 #flip i spin

		E_f = i * nn_sum

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
		# x, y = self.get_i()
		# x,y,_,_ = self.get_i_and_j()
		x,y = self.get_x_y()
		E_i, E_f = self.get_flip_E(x,y)
		if self.get_probability(-E_f+E_i) == True:
			# self.flip_spin(self.x,self.y)
			self.flip_spin(x,y)

	def Kawasaki(self):

		# self.get_i_and_j()
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

	def updatePlot(self,i):
		for i in range(self.N**2):
			self.update()
		self.plotFigure.clear()# Clear the old plot
		plt.imshow(self.array, interpolation = "nearest", cmap = "binary")# Make the new plot
		plt.axis('off')

	def Visualise(self):# Function that runs the animaion
		ani = animation.FuncAnimation(self.plotFigure, self.updatePlot)
		plt.show()


class Experiments():

	def __init__(self, N, dynamics):
		self.N = N
		self.Magnetisation = []
		self.Energy = []
		self.dynamics = dynamics
		self.init_state = "Uniform"

	def get_M(self,array):
		self.Magnetisation.append(np.sum(array))
		# return np.sum(array)
	def get_chi(self,M_vals,T):
		chi = 1/(self.N,T)
	def vary_T(self, min = 1, max = 3):
		# T_vals = np.linspace(min,max,(max-min)/0.1)
		T_vals = [2.7]

		for T in T_vals:

			sim = IsingModel(self.N, T, self.dynamics, self.init_state)

			for sweep in range(100): #equlibriate
				for i in range(self.N**2):
					sim.update()
				print "Equlibrium sweep", sweep
			#equlibriated

			for sweep in range(10000): #equlibriate
				for i in range(self.N**2):
					sim.update()
				if sweep%10 == 0:
					self.get_M(sim.array)
				print sweep
			#finished 1 T_val sample

		print self.Magnetisation
		print len(self.Magnetisation)
		print np.mean(self.Magnetisation)
				# else:
				# 	print "Pass"

def main():
	N = int(sys.argv[1])
	T = float(sys.argv[2])
	dynamics = sys.argv[3]
	init_state = sys.argv[4]

	# IsingModel(N,T,dynamics,init_state).Visualise()
	x = Experiments(N,dynamics)
	x.vary_T()
main()
