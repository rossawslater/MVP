import numpy as np

class Array(object):
    """Class to manage NxN Array"""
    def __init__(self, N):
        self.N = N
        self.make_array()

    def make_array(self):
        self.array = np.zeros((self.N,self.N))

    def find_NN(self,x,y):
        NN_array = np.zeros(9)
        NN_array[0] = (self.array[x - 1,y - 1])
        NN_array[1] =(self.array[x - 1, y])
        NN_array[2] =(self.array[x - 1,y + 1])
        NN_array[3] =(self.array[x ,y -1])
        NN_array[4] =(self.array[x ,y])
        NN_array[5] =(self.array[x ,y + 1])
        NN_array[6] =(self.array[x + 1,y-1])
        NN_array[7] =(self.array[x + 1,y])
        NN_array[8] =(self.array[x + 1,y + 1])
        NN_array.shape = (3,3)
        print NN_array
        return NN_array

class GoL():
    """Game of Life simulation"""
    def __init__(self, N):
        super().__init__(N)
        self.arg = arg

def main():
    array = Array(50)
    array.find_NN(5,5)

main()
