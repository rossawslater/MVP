import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from NbyNarray import create_array

class cahn_hilliard():

	def __init__(self, N, sweeps, dx, dt, alpha, kappa, mobility, phi_naught):

		self.sweep_number = sweeps
		self.N = N
		self.phi_naught = phi_naught
		self.phi_lattice = phi_naught*create_array.ones_array(N) + create_array.random_uniform_array(N, 1)/100

		self.dx = dx
		self.dt = dt
		self.alpha = alpha
		self.kappa = kappa
		self.mobility = mobility

		self.plot = plt.figure()


	def del_squared_phi(self, lattice):
		#Calculate array of laplacian values at each point.
		laplacian_array = (np.roll(lattice, 1, 1) + np.roll(lattice, -1, 1) + np.roll(lattice, 1, 0) + np.roll(lattice, -1, 0) - 4*lattice)/(self.dx**2)

		return(laplacian_array)


	def mu(self, lattice):
		#Calculate array of mu values at each point.
		mu_lattice = (-self.alpha * lattice) + (self.alpha*(lattice**3)) - self.kappa*self.del_squared_phi(lattice)

		return(mu_lattice)


	def update_phi(self, lattice):
		#updates all values of phi.
		mu_lattice = self.mu(lattice)
		increment_lattice = (self.mobility * self.dt) / (self.dx ** 2)  * (np.roll(mu_lattice, 1, 1) + np.roll(mu_lattice, -1, 1) + np.roll(mu_lattice, 1, 0) + np.roll(mu_lattice, -1, 0) - 4*mu_lattice)/(self.dx**2)
		future_lattice = lattice + increment_lattice

		return(future_lattice)


	def sweep(self):
		#Carry out update of full phi array.
		self.phi_lattice = create_array.set_boundaries(self.update_phi(self.phi_lattice), self.N)


	def sim(self):
		for i in range(self.sweep_number):
			self.sweep()


	def update_plot(self, i):
		self.plot.clear()
		self.sweep()
		plt.imshow(self.phi_lattice[1:self.N+1,1:self.N+1], cmap="gnuplot", origin = "lower", vmin=-1, vmax=1)


	def animate(self):
		ani = animation.FuncAnimation(self.plot, self.update_plot, interval=1)
		plt.show()


	def del_phi(self, lattice): # pass self.phi_lattice

		term1 = ((1/(2*self.dx))*(np.roll(lattice, 1, 0) - np.roll(lattice, -1, 0)))**2
		term2 = ((1/(2*self.dx))*(np.roll(lattice, 1, 1) - np.roll(lattice, -1, 1)))**2

		del_phi_lattice = (term1 + term2)

		return(del_phi_lattice)


	def get_free_energy(self, lattice): #pass self.phi_lattice

		free_energy_lattice = ((-self.alpha/2)*(lattice)**2) + ((self.alpha/4)*(lattice)**4) + ((self.kappa/2)*(self.del_phi(lattice)))
		total_free_energy = np.sum(free_energy_lattice)

		return(total_free_energy)


	def free_energy(self):
		n_iter = 20000

		data = np.zeros([n_iter, 2])

		for i in range(n_iter):
			data[i, 0] = i
			data[i, 1] = self.get_free_energy(self.phi_lattice)

			self.sweep()
			print(i)

		np.save('free_energy', data)


	def plot_free_energy(self):
		data = np.load('free_energy.npy')
		#print(data[:,0])
		plt.plot(data[:,0], data[:,1])
		plt.title("Free Energy over Time with $\phi$ = " + str(phi_naught))
		plt.xlabel("Number of Steps")
		plt.ylabel("Free Energy")
		plt.savefig('phinaughtZeroPointFive')

N = 100
sweeps = 25000
dx = 1
dt = 2.5
alpha = 0.1
mobility = 0.1
kappa = 0.1
phi_naught = 0.5

A = cahn_hilliard(N, sweeps, dx, dt, alpha, kappa, mobility, phi_naught)
A.animate()
#A.free_energy()
#A.plot_free_energy()
