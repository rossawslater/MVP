import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
# from scipy import ndimage
# import sys


class Array(object):
    """Class to manage the NxN Array"""
    def __init__(self, N):
        self.N = int(N)
        self.make_array()

    def make_array(self):
        self.array = np.zeros((self.N,self.N))
        return self.array

    def get_NN_array(self,x,y):
        NN_array = np.roll(np.roll(self.array,shift=-x+1,axis=0),shift=-y+1,axis=1) #returns 3x3 matrix centred on x,y
        return NN_array[:3,:3]

    def get_square_NNs(self,x,y):
        return [self.array[(x + 1)%self.N, y],self.array[x, (y + 1)%self.N],self.array[x-1, y],self.array[x, y-1]/
        self.array[(x + 1)%self.N, (y + 1)%self.N],self.array[x-1, (y + 1)%self.N],self.array[x-1, y-1],self.array[(x + 1)%self.N, y-1]]

    def get_plus_NNs(self,x,y):
        return [self.array[(x + 1)%self.N, y],self.array[x, (y + 1)%self.N],self.array[x-1, y],self.array[x, y-1]]

    def get_x_y(self):
		return np.random.randint(0,self.N), np.random.randint(0,self.N)

class GoL(Array):
    """Game of Life simulation"""
    def __init__(self, N):
        Array.__init__(self,N)

    def initalise_random(self):
        self.array = np.random.choice([0,1], (self.N,self.N))

    def get_NN(self,x,y):
        NN_array = self.get_NN_array(x,y)
        return np.sum(NN_array) - NN_array[1,1] #sums all Nearest Neighbours of central element

    def future_state(self,i,j,NNs): #just call current_state in funcion?
        current_state = self.array[i,j]
        if current_state == 1:
            if NNs <2:
                self.future_array[i,j] = 0
            if 2<= NNs <= 3:
                self.future_array[i,j] = 1
            if NNs > 3:
                self.future_array[i,j] = 0
        elif current_state == 0 and NNs == 3:
            self.future_array[i,j] = 1

    def update(self):
        self.future_array = np.copy(self.array)
        for i in range(self.N):
            for j in range(self.N):
                self.future_state(i,j,self.get_NN(i,j))
        self.array = np.copy(self.future_array)
        return self.array

    bee_hive = np.array([[0,1,0],[1,0,1],[1,0,1],[0,1,0]])

    oscillator = np.array([[1],[1],[1]])

    glider = np.array([[0,1,0],[0,0,1],[1,1,1]])

    def insert(self, state, x,y):
        temp = np.roll(np.roll(self.array,shift=-x,axis=0),shift=-y,axis=1)
        temp[0:len(state),0:len(state[0])] = state
        temp = np.roll(np.roll(temp,shift=x,axis=0),shift=y,axis=1)
        self.array = np.copy(temp)

class SIRS(Array):
    #0 is Susceptible, 1 is infected, 2 is recovered, 3 is immune
    def __init__(self, N, p1, p2, p3, immunity = False):
        Array.__init__(self,N)
        self.p1 = float(p1)
        self.p2 = float(p2)
        self.p3 = float(p3)
        self.immunity = float(immunity)
        self.initalise_random()

    def initalise_random(self):
        if self.immunity:
            # print "Immune"
            p = (1-self.immunity)/3
            self.array = np.random.choice([0,1,2,3], (self.N, self.N), p=[p,p,p,self.immunity])
        else:
            # print "no immunity"
            self.array = np.random.choice([0,1,2], (self.N, self.N))

    def initalise_uniform(self):
        self.array = np.zeros((self.N, self.N))
        self.array[self.N/2,self.N/2] = 1

    def update(self):
        for i in range(self.N**2):
            i,j = self.get_x_y()
            self.future_state(i,j)
        return self.array

    def future_state(self,i,j):

        current_state = self.array[i,j]
        r = np.random.rand()
        NN_array = self.get_plus_NNs(i,j)

        if current_state == 0:
            if 1 in NN_array:
                #if infected neighbor
                if r <= self.p1:
                    self.array[i,j] = 1

        elif current_state == 1:
            if r <= self.p2:
                self.array[i,j] = 2

        elif current_state == 2:
            if r <= self.p3:
                self.array[i,j] = 0

        elif current_state == 3:
            pass

class Visualise():
    def __init__(self, func):
        self.func = func
        self.plotFigure = plt.figure()
        self.show()

    def updatePlot(self,i): #function updates plot every 10 sweeps, as like measurements
        self.plotFigure.clear()# Clear the old plot
        plt.imshow(self.func(), interpolation = "nearest")#, cmap = "binary")# Make the new plot
        plt.axis('off')
        plt.colorbar()

    def show(self):# Function that runs the animaion
		ani = animation.FuncAnimation(self.plotFigure, self.updatePlot)
		plt.show()
