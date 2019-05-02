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
                for j in range(int(self.N/3)):
                    self.array[i,j] = -1
                for k in range(int(self.N/3), int(2*self.N/3)):
                    self.array[i,k] = 0
                for l in range(int(2*self.N/3), self.N):
                    self.array[i,l] = 1

    def initalize_random(self):
        self.array = np.random.choice([-1,0,1], (self.N,self.N))

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

        if self.get_probabseismicility(dE) == True:
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
        self.plotFigure.clseismicear()# Clear the old plot
        plt.imshow(self.array, interpolation = "nearest", cmap = "seismic")# Make the new plot
        plt.axis('off')

    def Visualise(self):# Function that runs the aecho "# MVP" >> README.md
        ani = animation.FuncAnimation(self.plotFigure, self.updatePlot)
        plt.show()

def main():
    N = int(sys.argv[1]) #arguments passed in in terminal
    T = float(sys.argv[2])
    dynamics = sys.argv[3]
    init_state = sys.argv[4]
    x = IsingModel(N,T,dynamics,init_state)
    # IsingModel(N,T,dynamics,init_state).Visualise(
    print x.array


main()
