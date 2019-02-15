import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

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

    def initalise_random(self):
        for x in range(0, self.N):
			for y in range(0, self.N):
				self.array[x,y] = np.random.choice([0,1])

class GoL(Array):
    """Game of Life simulation"""
    def __init__(self, N):
        Array.__init__(self,N)

    def get_NN(self,x,y):
        NN_array = self.get_NN_array(x,y)
        return np.sum(NN_array) - NN_array[1,1] #sums all Nearest Neighbours

    def future_state(self,i,j,current_state,NNs): #just call current_state in funcion?
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
                self.future_state(i,j,self.array[i,j],self.get_NN(i,j))
        self.array = np.copy(self.future_array)
        return self.array

class SIRS(Array):

    def __init__(self,N, p1, p2, p3):
        Array.__init__(self,N)
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3

    def update(self):
        self.future_array = np.copy(self.array)
        for i in range(self.N):
            for j in range(self.N):
                self.future_state(i,j,self.array[i,j],self.get_NN_array(i,j))
        self.array = np.copy(self.future_array)
        return self.array

    def future_state(self, i,j, current_state, NNs):
        r = np.random.rand()
        # if current_state == 1 and isin(1,self.array) == True:
        if 1 in NNs:
            print "infected Neighbour"
            if r < self.p1:
                self.future_array[i,j] = 0
        if current_state == 0:
            print ""
            if r < self.p2:
                self.future_array[i,j] = -1
        if current_state == -1:
            print ""
            if r < self.p3:
                self.future_array[i,j] = 1

class Visualise():
    def __init__(self, func):
        self.func = func
        self.plotFigure = plt.figure()
        self.show()

    def updatePlot(self,i): #function updates plot every 10 sweeps, as like measurements
		# self.func()
		self.plotFigure.clear()# Clear the old plot
		plt.imshow(self.func(), interpolation = "nearest")#, cmap = "binary")# Make the new plot
		plt.axis('off')

    def show(self):# Function that runs the animaion
		ani = animation.FuncAnimation(self.plotFigure, self.updatePlot)
		plt.show()

def main():
    # array = Array(50)
    # array.get_NN_array(5,5)
    # x = GoL(50)
    # x.initalise_random()
    y = SIRS(50,0.5,0.5,0.5)
    y.initalise_random()
    # print x.update()
    # x.Visualise()
    Visualise(y.update)
    # print x.get_NN_array(5,5)
    # print x.get_NN(5,5)

main()
