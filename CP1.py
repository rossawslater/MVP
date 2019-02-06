import numpy as np
import time
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys

class IsingModel():
	#sys argumets N, T , dynamics
	def __init__(self):#, N, T, random = True, J = 1):
		self.N = int(sys.argv[1])
		self.T = float(sys.argv[2])
		self.J = 1
		self.plotFigure = plt.figure()

		if sys.argv[3] == "Random":
			self.initalize_random()

		if sys.argv[3] == "Uniform":
			self.initalize_uniform()

	def initalize_uniform(self):
		self.array = np.ones([self.N, self.N])
		self.test_array = np.copy(self.array)

	def initalize_random(self):
		self.array = np.zeros([self.N, self.N])
		for x in range(0, self.N):
			for y in range(0, self.N):
				self.array[x,y] = np.random.choice([-1,1])
		self.test_array = np.copy(self.array)

	def get_i(self):
		i = np.random.randint(0,self.N,(2))
		self.x = i[0]; self.y = i[1]
		return self.x, self.y

	def flip_spin(self):#, x, y):
		self.array[self.x, self.y] *= -1

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

		i *= -1

		E_f = i * nn_sum

		return E_i, E_f

	def get_probability(self, E):
		r = np.random.rand()
		p = np.exp(-E/self.T)

		if E <= 0:
			self.flip_spin()
		elif r <= p:
			self.flip_spin()
		else:
			pass

	def Glauber(self):
		x, y = self.get_i()
		E_i, E_f = self.get_flip_E(x,y)
		self.get_probability(-E_f+E_i)


	def swap_spins(self):

		self.test_array = np.copy(self.array)

		r = np.random.randint(0,self.N,(4))

		self.i_x = r[0]; self.i_y = r[1]
		self.j_x = r[2]; self.j_y = r[3]


		# self.array[self.i_x, self.i_y], self.array[self.j_x, self.j_y] = self.array[self.j_x, self.j_y], self.array[self.i_x, self.i_y]
		i = self.array[self.i_x, self.i_y]
		j = self.array[self.j_x, self.j_y]

		self.test_array[self.i_x, self.i_y] = j
		self.test_array[self.j_x, self.j_y] = i

		self.i = i#self.test_array[self.j_x, self.j_y]
		self.j = j#self.test_array[self.i_x, self.i_y]

	def Kawasaki(self):
		self.get_i()
		self.swap_spins()
		if self.i != self.j:
			E_i = self.get_E(self.i_x, self.i_y)
			E_j = self.get_E(self.j_x, self.j_y)
			self.get_probability(E_j + E_j)
		else:
			pass

	def updatePlot(self,i):
		for i in range (2500):
			# self.Kawasaki()
			self.Glauber()
		self.plotFigure.clear()# Clear the old plot
		plt.imshow(self.array, interpolation = "nearest", cmap = "binary")# Make the new plot

	def Visualise(self):# Function that runs the animaion
		ani = animation.FuncAnimation(self.plotFigure, self.updatePlot)
		plt.show()



# IsingModel(50,1).Visualise()
x = IsingModel()#50,2.7)
# x.make_array()
# for i in range(100000):
# 	x.Glauber()
# plt.imshow(x.array, interpolation = "nearest", cmap = "binary")
# plt.show()
#
x.Visualise()
