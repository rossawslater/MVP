import numpy as np

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
