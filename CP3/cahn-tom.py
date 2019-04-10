"""
Checkpoint 3, Cahn-Hilliard Equation.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys

class Emulsion:

    def __init__(self, n, m, initialisation):

        self.phi_lattice = np.zeros([n, m])
        self.next_step = np.zeros([n, m])
        self.e_lattice = np.zeros([n, m])
        self.n = n
        self.m = m

        self.dx = 1
        self.dt = 2.5
        self.num  = 0

        self.alpha = 0.1
        self.kappa = 0.1
        self.M = 0.1
        self.phi_nought = 0.5

        self.plotWindow = plt.figure()
        self.sweep_number = 100


        if(initialisation == "random"):
            self.phi_lattice += self.phi_nought
            for i in range(self.n):
                for j in range(self.m):
                    self.phi_lattice[i, j] += np.random.random() / 100
                    
        elif(initialisation == "drop"):
            if(len(sys.argv) == 5):
                size = int(sys.argv[4])
            else:
                size = int(input("Define size (integer): "))
            
            for i in range(int(self.n / 2 - size / 2), int(self.n / 2 + size / 2)):
                for j in range(int(self.m / 2 - size / 2), int(self.m / 2 + size / 2)):
                    self.phi_lattice[i, j] = 1
            
        else:
            print("Invalid initialisation.")
            exit()


    def bounds(self, i, j):
        """
        returns position of (up, down, left, right) using boundary conditions
        """
        up = i + 1
        down = i - 1
        left = j - 1
        right = j + 1

        if(i == 0):
            down = self.n - 1
        if(i == self.n - 1):
            up = 0
        if(j == 0):
            left = self.m - 1
        if(j == self.m - 1):
            right = 0

        return up, down, left, right


    def del_squared_phi(self, i, j):
        """
        returns the laplacian of phi at X = (i, j)
        """
        up, down, left, right = self.bounds(i, j)

        factor = 1 / (self.dx ** 2)
        phi_sum_1 = self.phi_lattice[up, j] + self.phi_lattice[down, j]
        phi_sum_2 = self.phi_lattice[i, right] + self.phi_lattice[i, left]
        phi_sum = phi_sum_1 + phi_sum_2 - 4 * self.phi_lattice[i, j]

        return factor * phi_sum


    def mu(self, i, j):
        """
        returns the value of the chemical potential for X = (i, j) and n
        """
        phis = (-self.alpha * self.phi_lattice[i, j]) + (self.alpha * (self.phi_lattice[i, j] ** 3)) 
        del2_phi = -self.kappa * self.del_squared_phi(i, j)

        return phis + del2_phi
        

    def next_phi(self, i, j):
        """
        increments value of phi at X = (i, j) using euler algorithm
        """

        up, down, left, right = self.bounds(i, j)

        curr = self.phi_lattice[i, j]

        m_factor = (self.M * self.dt) / (self.dx ** 2)
        mus_1 = self.mu(up, j) + self.mu(down, j) # mus for i + and i -
        mus_2 = self.mu(i, left) + self.mu(i, right) # mus for j + and j -
        mus = mus_1 + mus_2 - (4 * self.mu(i, j)) 

        increment = m_factor * mus

        self.next_step[i, j] = curr + increment
        

    def update_lattice(self):
        """
        updates phi_lattice to next step
        """
        for i in range(self.n):
            for j in range(self.m):
                self.next_phi(i, j) # increments phi at lattice pos (i, j)

        self.phi_lattice = self.next_step.copy()
        self.num += 1

    def update_n(self):
        """
        updates the lattice sweep_number times
        """
        for i in range(self.sweep_number):
            self.update_lattice()
        sys.stdout.write("Step: %s \r" %(self.num))
        sys.stdout.flush()

#-----------------------------------------------------------------------------------------------

    def update_plot(self, i):
        """
        updates plot in plot window
        """
        self.update_n()
        
        self.plotWindow.clear()
        plt.imshow(self.phi_lattice, cmap="coolwarm", origin = "lower", vmin=-1, vmax=1)

        
    def animate(self):
        """
        animates phi_lattice
        """
        ani = animation.FuncAnimation(self.plotWindow, self.update_plot, interval=100)
        plt.show()
        

#-----------------------------------------------------------------------------------------------
        
    def del_phi(self, i, j): # ask about this function
        """
        returns magnitude of del of phi at i, j
        """
        up, down, left, right = self.bounds(i, j) 
        term1 = 0.5 * np.sqrt((self.phi_lattice[up, j] - self.phi_lattice[down, j]) **2)
        term2 = 0.5 * np.sqrt((self.phi_lattice[i, left] - self.phi_lattice[i, right]) **2)
        
        return np.sqrt((term1 ** 2) + (term2 ** 2))


    def free_energy(self, i, j):
        """
        returns the free energy for a point on the lattice
        """
        term1 = (-self.alpha / 2) * (self.phi_lattice[i, j] ** 2)
        term2 = (self.alpha / 4) * (self.phi_lattice[i, j] ** 4)
        term3 = (self.kappa / 2) * (self.del_phi(i, j) ** 2)
        
        return term1 + term2 + term3


    def tot_free_energy(self):
        """
        returns the total free energy of the lattice
        """
        f_tot = 0

        for i in range(self.n):
            for j in range(self.m):
                f_tot += self.free_energy(i, j)
  
        return f_tot
    

    def free_energy_time(self):
        """
        produces data for free energy against time
        """
        numsteps = 200

        data = np.zeros([numsteps, 2])

        for i in range(numsteps):
            data[i, 0] = i*100
            data[i, 1] = self.tot_free_energy()
            self.update_n()
        
        np.savetxt("free_energy_time.txt", data)


#----------------------------------------------------------------------------------------------

def plot_energies():
    """
    plots free energies against step number
    """
    data = np.loadtxt("free_energy_time.txt")
    
    ys = data[:, 1]
    xs = data[:, 0]

    plt.plot(xs, ys)
    plt.title("Free Energy over Time, phi_0 = 0.5")
    plt.xlabel("Number of Steps")
    plt.ylabel("Free Energy")
    plt.show()

def main():

    if(len(sys.argv) == 1):
       print("Required arguments: n, m")
       exit()
    
    n = int(sys.argv[1])
    m = int(sys.argv[2])
    initialisation = "random"
    if(len(sys.argv) > 3):
        initialisation = str(sys.argv[3])

    emulsion = Emulsion(n, m, initialisation)

    if(initialisation == "random"):
        if(len(sys.argv) > 4):
            if(sys.argv[4] == "animate"):
                emulsion.animate()
            elif(sys.argv[4] == "collect"):
                emulsion.free_energy_time()
            elif(sys.argv[4] =="plot"):
                plot_energies()
        else:
            print("Specify animate, collect, or plot")
            exit()
main()
