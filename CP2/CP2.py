import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
# from Array import *

class Array(object):
    """Class to manage the NxN Array"""
    def __init__(self, N):
        self.N = N
        self.make_array()

    def make_array(self):
        self.array = np.zeros((self.N,self.N))
        return self.array

    def get_NN_array(self,x,y):
        NN_array =np.roll(np.roll(self.array,shift=-x+1,axis=0),shift=-y+1,axis=1) #returns 3x3 matrix centred on x,y
        return NN_array[:3,:3]

class GoL(Array):
    """Game of Life simulation"""
    def __init__(self, N):
        Array.__init__(self,N)

    def initalise_random(self):
        for x in range(0, self.N):
			for y in range(0, self.N):
				self.array[x,y] = np.random.choice([0,1])

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

    def __init__(self,N, p1, p2, p3):
        Array.__init__(self,N)
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3

    def initalise_random(self):
        for x in range(0, self.N):
			for y in range(0, self.N):
				self.array[x,y] = np.random.choice([-1,0,1])

    def update(self):
        self.future_array = np.copy(self.array)
        for i in range(self.N):
            for j in range(self.N):
                self.future_state(i,j,self.get_NN_array(i,j))
        self.array = np.copy(self.future_array)
        return self.array

    def future_state(self, i,j, NNs):
        current_state = self.array[i,j]
        r = np.random.rand()
        if 0 in NNs:
            if r < self.p1:
                self.future_array[i,j] = 0
        if current_state == 0:
            if r < self.p2:
                self.future_array[i,j] = -1
        if current_state == -1:
            if r < self.p3:
                self.future_array[i,j] = 1

class Visualise():
    def __init__(self, func):
        self.func = func
        self.plotFigure = plt.figure()
        self.show()

    def updatePlot(self,i): #function updates plot every 10 sweeps, as like measurements
		self.plotFigure.clear()# Clear the old plot
		plt.imshow(self.func(), interpolation = "nearest")#, cmap = "binary")# Make the new plot
		plt.axis('off')

    def show(self):# Function that runs the animaion
		ani = animation.FuncAnimation(self.plotFigure, self.updatePlot)
		plt.show()

def main():
    # x = GoL(50)
    # x.insert(x.bee_hive,5,25)
    # x.insert(x.oscillator,35,12)
    # x.insert(x.glider, 25, 25)
    # Visualise(x.update)

    y = SIRS(50,0.75,0.75,0.75)
    # y.initalise_random()
    y.array = np.ones((y.N,y.N))
    y.array[24,24] = 0

    Visualise(y.update)

main()
