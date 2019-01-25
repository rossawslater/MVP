import numpy as np 
import time
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class IsingModel():

	def __init__(self, N, T):
		self.N = N
		self.T = T
		self.plotFigure = plt.figure()
		self.make_array()
	
	def make_array(self): #generate NxN matrix of random 1s and 0s
		self.array = np.zeros([self.N, self.N])
		for x in range(0, self.N):
			for y in range(0, self.N):
				self.array[x,y] = np.random.choice([-1,1])
		self.test_array = np.copy(self.array)
	
	def get_i(self):
		i = np.random.randint(0,self.N,(2))
		self.x = i[0]; self.y = i[1]
		
	def flip_spin(self, x, y):
		self.test_array[self.x, self.y] = - self.array[self.x, self.y]

	def get_E(self, x, y):
		i = self.test_array[x, y]
		
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

		self.E = (self.test_array[x_plus_1, y]*i + self.test_array[x_minus_1,y]*i 
				+ self.test_array[x, y_plus_1]*i + self.test_array[x,y_minus_1]*i)

		return self.E


	def get_probability(self, E):
		if E <= 0:
			self.array = np.copy(self.test_array)
		if E > 0:
			r = np.random.rand()
			p = np.exp(-E/self.T)
			print r,p
			if r <= p:
				self.array = np.copy(self.test_array)
			else:
				pass


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

	def Glauber(self):
		self.get_i()
		E_1 = self.get_E(self.x, self.y)
		self.flip_spin(self.x, self.y)
		E_2 = self.get_E(self.x, self.y)
		self.get_probability(E_2-E_1)

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
	    self.Kawasaki()
	    self.plotFigure.clear()# Clear the old plot
	    plt.imshow(self.array, interpolation = "nearest", cmap = "binary")# Make the new plot 
	
	def Visualise(self):# Function that runs the animaion 
		# self.make_array()
		ani = animation.FuncAnimation(self.plotFigure, self.updatePlot)#, interval=1 )
		plt.show()


	
# IsingModel(50,1).Visualise()
x = IsingModel(50,0.01)
# x.make_array()
# for i in range(1000000):
# 	x.Kawasaki()
# plt.imshow(x.array, interpolation = "nearest", cmap = "binary")
# plt.show()

x.Visualise()