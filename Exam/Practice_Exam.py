import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random
import sys

class CahnHillard():

    def __init__(self, N):
        self.N = int(N)
        # self.phi = np.random.normal(float(phi0), 0.1, (self.N, self.N))
        self.phi = np.random.uniform(-0.1,0.1,(self.N,self.N)) + 0.5

        self.sigma = 10
        self.k = 0.01
        self.D = 1.0
        self.rho = np.zeros((self.N, self.N))
        self.r_vals = np.zeros((self.N, self.N))


        self.dt = 0.01
        self.dx = 1.0

        self.phi_list= []

    def get_rho_array(self):
        for i in range(self.N):
            for j in range(self.N):
                r = np.sqrt((i - self.N/2)**2 + (j - self.N/2)**2)
                self.r_vals[i,j] = r
                rho = np.exp(-r**2 / self.sigma**2)
                self.rho[i,j] = rho

    def get_phi_array(self):
        self.new_phi_array = self.phi + (self.D*self.dt/(self.dx**2)) * (np.roll(self.phi, 1,0) + np.roll(self.phi, -1,0) + np.roll(self.phi, 1,1) + np.roll(self.phi, -1,1)\
        - 4.0 * self.phi) + self.rho - self.k * self.phi
        self.phi = np.copy(self.new_phi_array)
        return self.phi

    def update(self):
        # for sweep in range(5000):
            # self.get_mu_array()
        self.get_phi_array()
        self.phi_list.append(np.mean(self.phi))
        return self.phi

    def get_f(self):
        #----------Using Central Difference Method
        mx = (np.roll(self.Du_array, 1,0) - np.roll(self.Du_array, -1,0))/(2*self.dx)
        my = (np.roll(self.Du_array, 1,1) - np.roll(self.Du_array, -1,1))/(2*self.dx)
        grad = np.sqrt(mx**2 + my**2)
        self.f = -self.a/2.0*(self.phi**2) + (self.a/4.0)*(self.phi**4) + (self.k/(2)*(grad**2))
        self.all_f.append(np.sum(self.f))
        return np.sum(self.f)

    def plot_f(self):
        plt.plot(self.all_f)
        plt.title("Free Energy against Time")
        plt.xlabel("100 Sweeps")
        plt.ylabel("Free Energy [Arbritary Units]")
        plt.show()

class Visualise():
    def __init__(self, func):
        self.func = func
        self.plotFigure = plt.figure()
        self.show()

    def updatePlot(self,i): #function updates plot every 10 sweeps, as like measurements
        self.plotFigure.clear()# Clear the old plot
        plt.imshow(self.func(), cmap = "seismic")#, interpolation = "nearest")#, cmap = "binary")# Make the new plot
        plt.axis('off')
        plt.colorbar()

    def show(self):# Function that runs the animaion
        ani = animation.FuncAnimation(self.plotFigure, self.updatePlot)
        plt.show()

def main():
    N = sys.argv[1]
    # phi0 = sys.argv[2]

    x = CahnHillard(N)
    x.get_rho_array()
    # print x.rho
    # print x.phi
    # plt.imshow(x.rho)
    # plt.show()
    for i in range(1000):
        x.update()
    plt.plot(x.phi_list)
    plt.show()
    plt.scatter(x.r_vals, x.phi)
    plt.show()
    # Visualise(x.update)
main()
