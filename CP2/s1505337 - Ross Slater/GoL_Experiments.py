from CP2 import *
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage

class Experiments():

    def __init__(self, N = 50):
        self.N = N

    def get_CoM(self, array):
        return ndimage.measurements.center_of_mass(array)

    def glider_com(self):
        self.sim = GoL(self.N)
        self.sim.insert(self.sim.glider,0,0)
        COMs = []
        COM_x = []
        COM_y = []
        sim_len = 100
        for i in range(sim_len):
            self.sim.update()
            COMs.append(self.get_CoM(self.sim.array))
            x, y = self.get_CoM(self.sim.array)
            COM_x.append(x)
            COM_y.append(y)
        COM_x = np.array(COM_x)
        COM_y = np.array(COM_y)
        pos = np.sqrt(COM_x**2 + COM_y**2)
        np.savetxt(str("Data/COM_Pos_Glider.txt"), pos)
        speed = np.polynomial.polynomial.polyfit(range(100), pos, 1)[1]
        return COM_x, COM_y, pos, speed

def main():
    x, y, pos, speed = Experiments().glider_com()
    print speed
main()
