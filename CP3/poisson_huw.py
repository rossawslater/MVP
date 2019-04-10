import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

class poisson():

    def __init__(self, n_iter, N, error, initialisation):

        self.omega = 1.75
        self.n_iter = n_iter
        self.error = error
        self.N = N
        self.mid = int(N/2)

        self.ro_lattice = np.zeros([self.N+2, self.N+2, self.N+2])


        if initialisation == 'point' or 'point charge' or 'Point' or 'Point Charge':
            self.ro_lattice[self.mid+2, self.mid+2, self.mid+2] = 1
            print('Point charge initialised.')

        elif initialisation == 'Wire' or 'wire':
            self.ro_lattice[self.mid, self.mid, :] = 1
            print('Wire initialised.')

        self.phi_lattice = np.zeros([self.N+2, self.N+2, self.N+2])


    def jacobi(self, lattice):
        #pass phi lattice, return jacobi lattice

        right = np.roll(lattice, 1, 0)
        left = np.roll(lattice, -1, 0)
        up = np.roll(lattice, 1, 1)
        down = np.roll(lattice, -1, 1)
        back = np.roll(lattice, 1, 2)
        forward = np.roll(lattice, -1, 2)

        return((1 / 6)*(up + down + left + right + forward + back + self.ro_lattice))


    def gauss_seidel(self, lattice):

        update_lattice = np.copy(lattice)

        for i in range(1,self.N+1):
            for j in range(1,self.N+1):
                for k in range(1,self.N+1):
                    gauss_seidel = 1/6.0*(lattice[(i+1),j,k]+update_lattice[(i-1),j,k]+lattice[i,(j+1),k]+update_lattice[i,(j-1),k]+lattice[i,j,(k+1)]+update_lattice[i,j,(k-1)]+self.ro_lattice[i,j,k])

                    update_lattice[i,j,k] = (1-self.omega)*lattice[i,j,k] + self.omega*gauss_seidel

        return(update_lattice)


    def original_gauss_seidel(self, lattice):

        update_lattice = np.copy(lattice)

        for i in range(self.N):
            for j in range(self.N):
                for k in range(self.N):
                    update_lattice[i,j,k] = 1/6*(lattice[(i+1)%self.N,j,k]+update_lattice[(i-1)%self.N,j,k]+lattice[i,(j+1)%self.N,k]+update_lattice[i,(j-1)%self.N,k]+lattice[i,j,(k+1)%self.N]+update_lattice[i,j,(k-1)%self.N]+self.ro_lattice[i,j,k])

        return(update_lattice)


    def E_field(self, lattice): # pass self.phi_lattice

        Ex_field_lattice, Ey_field_lattice, Ez_field_lattice = np.gradient(lattice[:, :, :])

        return(Ex_field_lattice, Ey_field_lattice, Ez_field_lattice)


    def B_field(self, lattice):

        bx_field_lattice, by_field_lattice = np.gradient(lattice[:, :, self.mid])

        return(bx_field_lattice, by_field_lattice)


    def update_step(self, i):

        #update_lattice = self.jacobi(self.phi_lattice)
        update_lattice = self.gauss_seidel(self.phi_lattice)
        #update_lattice = self.original_gauss_seidel(self.phi_lattice)
        self.check_error(self.phi_lattice, update_lattice, i)
        self.phi_lattice = np.copy(update_lattice)
        self.bound_zero(self.phi_lattice)



    def bound_zero(self, lattice):
        #sets bounaries of lattice to 0

        lattice[0, :, :] = 0
        lattice[-1, :, :] = 0
        lattice[:, 0, :] = 0
        lattice[:, -1, :] = 0
        lattice[:, :, 0] = 0
        lattice[:, :, -1] = 0

        return(lattice)


    def run(self):

        for i in range(self.n_iter):
            self.update_step(i)


    def check_error(self, phi_lattice, update_lattice, i):
        #determines whether next_step and phi_lattice are the same within an error, ends sim if so

        diff_lattice = (phi_lattice[1:self.N+1,1:self.N+1,1:self.N+1]- update_lattice[1:self.N+1,1:self.N+1,1:self.N+1])
        diff = np.max(abs(diff_lattice))
        #print(diff)
        if(diff <= self.error):
            print(i)
            E = self.E_field(phi_lattice)
            B = self.B_field(phi_lattice)
            np.save('E_x', E[0])
            np.save('E_y', E[1])
            np.save('E_z', E[2])
            np.save('B_x', B[0])
            np.save('B_y', B[1])
            np.save('Electrostatic_Potential', phi_lattice)
            exit()


def plot_potential(lattice):


    norm = cm.colors.Normalize(vmax=np.max(lattice), vmin=0)
    plt.imshow(lattice[:, :, 50], interpolation = 'none', cmap = 'gnuplot', norm=norm)
    plt.title('Potential Strength Cut Through Mid-Plane')
    plt.xlabel('x-axis')
    plt.ylabel('y-axis')
    plt.colorbar()
    plt.show()

    plt.plot(np.arange(100), lattice[:, 50, 50])
    plt.title('Electrostatic Potential Along x-axis')
    plt.xlabel('x-axis')
    plt.ylabel('Electrostatic Potential')
    plt.savefig('Electrostatic Potential Along x-axis')

    plt.plot(np.arange(100), lattice[50, :, 50])
    plt.title('Electrostatic Potential Along y-axis')
    plt.xlabel('y-axis')
    plt.ylabel('Electrostatic Potential')
    plt.savefig('Electrostatic Potential Along y-axis')

    plt.plot(np.arange(100), lattice[50, 50, :])
    plt.title('Electrostatic Potential Along z-axis')
    plt.xlabel('z-axis')
    plt.ylabel('Electrostatic Potential')
    plt.savefig('Electrostatic Potential Along z-axis')

def plot_E():

    x_lattice = np.load('E_x.npy')
    y_lattice = np.load('E_y.npy')
    z_lattice = np.load('E_z.npy')

    E_field_mag_lattice = np.sqrt(x_lattice**2 + y_lattice**2 + z_lattice**2)

    axis = np.arange(100)

    plt.plot(np.arange(100), E_field_mag_lattice[:, 50, 50])
    plt.title('Electric Field Strength Along x-axis')
    plt.xlabel('x-axis')
    plt.ylabel('Electrostatic Field Strength')
    plt.savefig('Electric Field Strength Along x-axis')

    plt.plot(np.arange(100), E_field_mag_lattice[50, :, 50])
    plt.title('Electric Field Strength Along y-axis')
    plt.xlabel('y-axis')
    plt.ylabel('Electrostatic Field Strength')
    plt.savefig('Electric Field Strength Along y-axis')

    plt.plot(np.arange(100), E_field_mag_lattice[50, 50, :])
    plt.title('Electric Field Strength Along z-axis')
    plt.xlabel('z-axis')
    plt.ylabel('Electrostatic Field Strength')
    plt.savefig('Electric Field Strength Along z-axis')

def quiver_E():

    x = np.arange(100)
    y = np.arange(100)
    z = np.arange(100)

    u_array = np.load('E_x.npy')
    v_array = np.load('E_y.npy')
    w_array = np.load('E_z.npy')

    u = u_array[:, :, 50]
    v = v_array[:, :, 50]
    w = w_array[:, :, 50]

    fig, ax = plt.subplots()
    ax.quiver(x, y, u, v, pivot='mid')

    plt.show()

def quiver_B():

    x = np.arange(100)
    y = np.arange(100)

    u_array = np.load('B_x.npy')
    v_array = np.load('B_y.npy')

    u = u_array
    v = -v_array

    fig, ax = plt.subplots()
    ax.quiver(x, y, u, v, pivot='tip')

    plt.show()


def main():

    n_iter = 10000
    N = 50
    error = 0.00001
    initialisation = 'Point Charge'

    A = poisson(n_iter, N, error, initialisation)
    A.run()


# main()
# plot_E()
# quiver_E()
# quiver_B()
#plot_potential(np.load('Electrostatic_Potential.npy'))
