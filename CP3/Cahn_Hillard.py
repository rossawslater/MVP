import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random
import sys

class CahnHillard():

    def __init__(self, N, phi0):
        self.N = int(N)
        self.phi = np.random.normal(float(phi0), 0.0001, (self.N, self.N))
        # self.phi = np.random.choice([-0.5,0.5],size = (self.N, self.N))
        self.a = 0.1
        self.b = 0.1
        self.k = 0.1
        self.M = 0.1

        self.dt = 0.75
        self.dx = 0.75

        self.all_f = []

    def get_mu_array(self):
        self.mu_array = (-self.a * self.phi) + (self.b * self.phi**3) - (self.k/(self.dx**2) * (np.roll(self.phi, 1,0) + np.roll(self.phi, -1,0) + np.roll(self.phi, 1,1) + np.roll(self.phi, -1,1) - 4.0*self.phi))

    def get_phi_array(self):
        self.new_phi_array = self.phi + (self.M*self.dt/(self.dx**2) * (np.roll(self.mu_array, 1,0) + np.roll(self.mu_array, -1,0) + np.roll(self.mu_array, 1,1) + np.roll(self.mu_array, -1,1) - 4.0*self.mu_array))
        self.phi = np.copy(self.new_phi_array)
        return self.phi

    def update(self):
        for sweep in range(1000):
            self.get_mu_array()
            self.get_phi_array()
        print self.get_f()
        return self.phi

    def get_f(self):
        #----------Using Central Difference Method
        mx = (np.roll(self.mu_array, 1,0) - np.roll(self.mu_array, -1,0))/(2*self.dx)
        my = (np.roll(self.mu_array, 1,1) - np.roll(self.mu_array, -1,1))/(2*self.dx)
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
    phi0 = sys.argv[2]
    Visualise(CahnHillard(N, phi0).update)
main()
