"""
Checkpoint3, Poisson's equation
"""
import numpy as np
import sys 
import matplotlib.pyplot as plt

class Poisson:

    def __init__(self, n, m, o, initialisation="point", margin=0.01):
        
        self.n = n
        self.m = m
        self.o = o
        
        self.ro_lattice = np.zeros([n, m, o]) # 'input lattice'

        self.phi_lattice = np.zeros([n, m, o])
        self.next_step = np.zeros([n, m, o])

        numpoints = m * n * o
        self.e_field = np.zeros([numpoints, 6]) # e stored as x y z xdir ydir zdir

        self.margin = margin

        if(initialisation == "point"):
            i = int(self.n / 2)
            j = int(self.m / 2)
            k = int(self.o / 2)
            self.ro_lattice[i, j, k] = 1


    def jacobi(self, i, j, k): 
        """
        returns jacobi update of point i, j, k
        """
        i_terms = self.phi_lattice[i + 1, j, k] + self.phi_lattice[i - 1, j, k]
        j_terms = self.phi_lattice[i, j + 1, k] + self.phi_lattice[i, j - 1, k]
        k_terms = self.phi_lattice[i, j, k + 1] + self.phi_lattice[i, j, k + 1]
        ro_term = self.ro_lattice[i, j, k]

        return (1 / 6) * (i_terms + j_terms + k_terms + ro_term)


    def update_step(self): # calls jacobi
        """
        updates lattice to next step
        """
        self.next_step = np.zeros([self.n, self.m, self.o])
        
        for i in range(self.n - 1):
            for j in range(self.m - 1):
                for k in range(self.o - 1):
                    self.next_step[i, j, k] = self.jacobi(i, j, k) # calc next step
        
        self.end_iteration() # ends if same as last step

        self.phi_lattice = self.next_step.copy()
        self.bound_zero()


    def bound_zero(self):
        """
        sets bounaries of lattice to 0, dirichlet boundaries, x = 0 instead of dx=0
        """
        self.phi_lattice[:, :, 0] = 0
        self.phi_lattice[:, :, -1] = 0
        self.phi_lattice[0, :, :] = 0
        self.phi_lattice[-1, :, :] = 0
        self.phi_lattice[:, 0, :] = 0
        self.phi_lattice[:, -1, :] = 0


    def distance(self, pos1, pos2):
        """
        returns distance between 2 points as a vector, given arrays of [i, j, k]
        """
        rx = (pos2[0] - pos1[0]) ** 2
        ry = (pos2[1] - pos1[1]) ** 2
        rz = (pos2[2] - pos1[2]) ** 2
        
        return rx, ry, rz
        
#------------------------------------------------------------------------------------------------     
    def e_comp(self, pos1, pos2):
        """
        returns component of e at pos1 from pos2 as ex, ey, ez
        """
        rx, ry, rz = self.distance(pos1, pos2)
        ro = self.ro_lattice[pos2[0], pos2[1], pos2[2]]
        factor = ro / (4 * np.pi)

        if(rx == 0):
            ex = 0
        else:
            ex = factor / (rx ** 2)
        if(ry == 0):
            ey = 0
        else:
            ey = factor / (ry ** 2)
        if(rz == 0):
            ez = 0
        else:
            ez = factor / (rz ** 2)
        
        return ex, ey, ez


    def e_comps(self, x, y, z):
        """
        returns the components of the e field at x, y, z
        """
        e_pos = np.array([x, y, z])
        extot = 0
        eytot = 0
        eztot = 0

        for i in range(self.n):
            for j in range(self.m):
                for k in range(self.o):
                    ro_pos = np.array([i, j, k])
                    ex, ey, ez = self.e_comp(e_pos, ro_pos) # find e at e_pos from ro_pos
                    extot += ex
                    eytot += ey
                    eztot += ez
        
        return extot, eytot, eztot


    def e_comps_point(self, x, y, z):
        """
        returns comps of electric field from point charge at centre 
        """
        e_pos = np.array([x, y, z])

        i = int(self.n / 2)
        j = int(self.m / 2)
        k = int(self.o / 2)

        ro_pos = np.array([i, j, k])
        ex, ey, ez = self.e_comp(e_pos, ro_pos)

        return ex, ey, ez


    def find_e_field(self):
        """
        finds the e field at all points for timestep n
        """
        n = 0
        num = self.m*self.n*self.o
        
        for i in range(self.n):
            for j in range(self.m):
                for k in range(self.o):
                    self.e_field[n, 0] = i
                    self.e_field[n, 1] = j
                    self.e_field[n, 2] = k
                    ex, ey, ez = self.e_comps_point(i, j, k)
                    self.e_field[n, 3] = ex
                    self.e_field[n, 4] = ey
                    self.e_field[n, 5] = ez
                    n += 1
                    sys.stdout.write("Finding E-field: %.1f%% \r" %((n / num)*100))
                    sys.stdout.flush()

#-----------------------------------------------------------------------------------------------
    
    def collect_data(self):
        """
        overall data collection
        """
        steps = 0
        max_steps = 500
        while(steps < max_steps):
            self.update_step()
            steps += 1
            if(steps%100 == 0):
                sys.stdout.write("%s \r" %(steps))
                sys.stdout.flush()
        self.end()

#--------------------------------------------------------------------------------------------

    def end_iteration(self):
        """
        determines whether next_step and phi_lattice are the same within an error, ends sim if so
        """
        diff_lattice = np.sqrt((self.phi_lattice - self.next_step) ** 2)
        tot_diff = np.sum(diff_lattice)
        if(tot_diff <= self.margin):
            self.end()


    def end(self):
        """
        ends iterations and returns values
        """
        sys.stdout.write("\n")
        numpoints = self.m * self.n * self.o
        phi_save = np.zeros([numpoints, 4])
        n = 0
        for i in range(self.n):
            for j in range(self.m):
                for k in range(self.o):                    
                    phi_save[n, 0] = i
                    phi_save[n, 1] = j
                    phi_save[n, 2] = k
                    phi_save[n, 3] = self.phi_lattice[i, j, k]
                    n += 1
       
        np.savetxt("phi_field.txt", phi_save)
        
        self.find_e_field()
        np.savetxt("e_field.txt", self.e_field)
        exit()

#-------------------------------------------------------------------------------------------

def get_axis(axis):
    """
    takes string axis name, returns int and limit
    """
    if(axis == "x"):
        ax = 0
        lim = None
    elif(axis =="y"):
        ax = 1
        lim = 899 # after lim, repeats same values
    elif(axis == "z"):
        ax = 2
        lim = 29
    else:
        print("invalid axis")
        exit()

    return ax, lim


def r_e_plot():
    """
    plots 3 dimension with e components
    """
    array = np.loadtxt("e_field.txt")

    x1 = array[:, 0] # x direction
    y1 = array[:, 3]

    ylim = sys.argv[3] * sys.argv[4]
    x2 = array[:ylim, 1] # y direction
    y2 = array[:ylim, 4]
    
    zlim = sys.argv[4] - 1
    x3 = array[:zlim, 2] # z direction
    x4 = array[:zlim, 5]

    plt.subplot(3,1,1)
    plt.plot(x1, y1)
    plt.title("Change in Ex with x")
    plt.ylabel("E-field in x direction")
    plt.xlabel("x")

    plt.subplot(3,1,2)
    plt.plot(x2, y2)
    plt.title("Change in Ey with y")
    plt.ylabel("E-field in y direction")
    plt.xlabel("y")

    plt.subplot(3,1,3)
    plt.plot(x3, y3)
    plt.title("Change in Ez with z")
    plt.ylabel("E-field in z direction")
    plt.xlabel("z")

    plt.show()


def r_phi_plot(axis):
    """
    plots single dimension with phi
    """
    ax, lim = get_axis(axis)

    array = np.loadtxt("phi_field.txt")

    xs = array[:lim, ax]
    ys = array[:lim, 3]

    plt.plot(xs, ys)
    plt.title("Change in phi with " + axis)
    plt.ylabel("Phi")
    plt.xlabel(axis)
    plt.show()


def main():
    
    method = sys.argv[1]
    n = int(sys.argv[2])
    m = int(sys.argv[3])
    o = int(sys.argv[4])

    if(method == "collect"):
        if(len(sys.argv) > 6):
            initialisation = sys.argv[5] 
            margin = float(sys.argv[6])
            p = Poisson(n, m, o, initialisation, margin)
            p.collect_data()
        else:
            print("Initialisation and margin of error required.")

    if(method == "plot_e"):
        r_e_plot()

    if(method == "plot_phi"):
        if(len(sys.argv) > 5):
            axis = sys.argv[5]
            r_phi_plot(axis)

    if(method == "test"):
        initialisation = sys.argv[5] 
        margin = float(sys.argv[6])
        p = Poisson(n, m, o, margin)
        p.end_iteration()

main()


#PROBLEMS: E-FIELD=0 EVERYWHER
