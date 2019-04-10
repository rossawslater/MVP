import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random
import sys

class Poisson():

    def __init__(self, N, Accuracy = 0.00001):
        self.N = int(N) + 2
        self.dx = 1
        self.rho = np.zeros((self.N, self.N, self.N))
        # self.rho = np.random.uniform(0.1, -0.1, size = (self.N, self.N, self.N))
        self.phi = np.zeros((self.N, self.N, self.N))
        self.tolerance = float(Accuracy)
        self.complete = False
        self.omega = 1

    def place_central_charge(self):
        mid = int(len(self.rho[0])/2)
        self.rho[mid,mid,mid] = 1

    def place_wire(self):
        mid = int(len(self.rho[0])/2)
        self.rho[mid,mid,:] = 1

    def get_new_phi_array(self):
        self.old_phi = np.copy(self.phi)

        self.new_phi_array = (np.roll(self.phi, 1,0) + np.roll(self.phi, -1,0) \
        + np.roll(self.phi, 1,1) + np.roll(self.phi, -1,1) \
        + np.roll(self.phi, 1,2) + np.roll(self.phi, -1,2))
        self.new_phi_array += self.dx**2 * self.rho #still keeping dx^2 term even though set to one earlier
        self.new_phi_array /= 6

        self.set_boundaries()
        self.phi = np.copy(self.new_phi_array)
        return self.phi

    def set_boundaries(self):
        self.new_phi_array[0,:,:] = 0
        self.new_phi_array[-1,:,:] = 0
        self.new_phi_array[:,0,:] = 0
        self.new_phi_array[:,-1,:] = 0
        self.new_phi_array[:,:,0] = 0
        self.new_phi_array[:,:,-1] = 0

    def iterate_phi(self, i, j, k):
        new_phi = (self.old_phi[i+1,j,k] + self.old_phi[i-1,j,k] + self.old_phi[i,j+1,k] + self.old_phi[i,j-1,k] + self.old_phi[i,j,k+1] + self.old_phi[i,j,k-1] + self.rho[i,j,k])/6
        self.phi[i,j,k] = new_phi

    def iterate_phi_GS(self, i, j, k):
        new_phi = (self.phi[i+1,j,k] + self.phi[i-1,j,k] + self.phi[i,j+1,k] + self.phi[i,j-1,k] + self.phi[i,j,k+1] + self.phi[i,j,k-1] + self.rho[i,j,k])/6
        self.phi[i,j,k] = (1-self.omega)*self.phi[i,j,k] + self.omega*new_phi

    def check_accuracy(self):
        if np.max(abs(self.phi - self.old_phi)) <= self.tolerance:
            print (np.max(abs(self.phi - self.old_phi)))
            self.complete = True

    def update_Jacobi(self):
        self.old_phi = np.copy(self.phi)
        self.get_new_phi_array()
        return self.phi

    def update_Jacobi_Iterative(self):
        self.old_phi = np.copy(self.phi)
        for i in range(1,self.N-1):
            for j in range(1,self.N-1):
                for k in range(1,self.N-1):
                    self.iterate_phi(i,j,k)

    def update_Gauss_Seidel(self):
        self.old_phi = np.copy(self.phi)
        for i in range(1,self.N-1):
            for j in range(1,self.N-1):
                for k in range(1,self.N-1):
                    self.iterate_phi_GS(i,j,k)

    def get_E(self):
        self.Ex, self.Ey, self.Ez = np.gradient(self.phi)
        np.savetxt("Ex_XY_wire.txt",self.Ex[:, int(self.N/2), :])
        np.savetxt("Ey_XY_wire.txt", self.Ey[:, int(self.N/2), :])
        # np.savetxt("EZ.txt", self.Ez)
        return self.Ex,self.Ey,self.Ez

    def get_B(self):
        self.Bx, self.By = np.gradient(self.phi[:,:,int(self.N/2)])

        np.savetxt("Bx_wire.txt",self.Bx)
        np.savetxt("By_wire.txt", self.By)
        return self.Bx, self.By

    def norm_fields(self):
        E_norm = np.sqrt(self.Ex**2 + self.Ey**2 + self.Ez**2)
        B_norm = np.sqrt(self.Bx**2 + self.By**2)
        self.Ex /= E_norm
        self.Ey /= E_norm
        self.Ez /= E_norm
        self.Bx /= B_norm
        self.By /= B_norm

#-------Plotting----------------------------------------------------------------

    def plot_electrostatic_potential_x(self):
        plt.plot(self.phi[:,int(self.N/2), int(self.N/2)])
        plt.title("Electrostatic Potential Through Midpoint Across x-axis")
        plt.xlabel("x-axis")
        plt.ylabel("Electrostatic Potential")
        plt.show()

    def plot_electrostatic_potential_y(self):
        plt.plot(self.phi[int(self.N/2),:,int(self.N/2)])
        plt.title("Electrostatic Potential Through Midpoint Across y-axis")
        plt.xlabel("y-axis")
        plt.ylabel("Electrostatic Potential")
        plt.show()

    def plot_electrostatic_potential_z(self):
        plt.plot(self.phi[int(self.N/2),int(self.N/2),:])
        plt.title("Electrostatic Potential Through Midpoint Across z-axis")
        plt.xlabel("z-axis")
        plt.ylabel("Electrostatic Potential")
        plt.show()

    def plot_electrostatic_potentials(self):
        self.plot_electrostatic_potential_x()
        self.plot_electrostatic_potential_y()
        self.plot_electrostatic_potential_z()
#------------------------------------------------------------------------------
    def plot_electric_field_x(self):
        plt.plot(abs(self.Ex[:,int(self.N/2), int(self.N/2)]))
        plt.title("Electric Field Through Midpoint Across x-axis")
        plt.xlabel("x-axis")
        plt.ylabel("Electric Field Strength")
        plt.show()

    def plot_electric_field_y(self):
        plt.plot(abs(self.Ey[int(self.N/2),:, int(self.N/2)]))
        plt.title("Electric Field Through Midpoint Across y-axis")
        plt.xlabel("y-axis")
        plt.ylabel("Electric Field Strength")
        plt.show()

    def plot_electric_field_z(self):
        plt.plot(abs(self.Ez[int(self.N/2), int(self.N/2),:]))
        plt.title("Electric Field Through Midpoint Across z-axis")
        plt.xlabel("z-axis")
        plt.ylabel("Electric Field Strength")
        plt.show()

    def plot_electric_fields(self):
        self.plot_electric_field_x()
        self.plot_electric_field_y()
        self.plot_electric_field_z()
#------------------------------------------------------------------------------

    def E_vector_plot_XY(self):

        U = (self.Ex[:, :, int(self.N/2)])
        V = (self.Ey[:, :, int(self.N/2)])

        fig, ax = plt.subplots()
        q = ax.quiver(V,U, pivot='mid', angles='xy')

        plt.title("Electric Field")
        plt.xlabel("x-axis")
        plt.ylabel("y-axis")
        plt.show()

        plt.title("Electric Field")
        plt.xlabel("x-axis")
        plt.ylabel("y-axis")
        plt.imshow(abs(self.phi[:,:,int(self.N/2)]), origin = "lower")
        plt.colorbar()

        plt.show()

    def E_vector_plot_XZ(self):

        U = (self.Ex[:, int(self.N/2), :])
        V = (self.Ez[:, int(self.N/2), :])

        fig, ax = plt.subplots()
        q = ax.quiver(V,U, pivot='mid', angles='xy')
        plt.title("Electric Field")
        plt.xlabel("x-axis")
        plt.ylabel("Z-axis")

        plt.imshow(abs(self.phi[:,int(self.N/2),:]), origin = "lower")
        plt.colorbar()

        plt.show()

    def E_vector_plot_YZ(self):

        U = (self.Ey[int(self.N/2),:,:])
        V = (self.Ez[int(self.N/2),:,:])

        fig, ax = plt.subplots()
        q = ax.quiver(V,U, pivot='mid', angles='xy')
        plt.title("Electric Field")
        plt.xlabel("y-axis")
        plt.ylabel("z-axis")
        plt.imshow(abs(self.phi[int(self.N/2),:,:]), origin = "lower")
        plt.colorbar()

        plt.show()
#------------------------------------------------------------------------------
    def B_field_plot(self):
        plt.title("Magnetic Field Strength")
        plt.xlabel("x-axis")
        plt.ylabel("y-axis")
        plt.imshow(abs(self.phi[:,:,int(self.N/2)]), origin = "lower")
        plt.colorbar()
        plt.show()

    def B_vector_plot_XY(self):
        fig, ax = plt.subplots()
        q = ax.quiver(self.Bx,-self.By, pivot='tip', angles='xy')
        plt.title("Magnetic Field")
        plt.xlabel("x-axis")
        plt.ylabel("y-axis")
        # plt.imshow(abs(self.phi[:,:,int(self.N/2)]), origin = "lower")
        # plt.colorbar()

        plt.show()
#------------------------------------------------------------------------------
    def run(self):
        i = 0
        while self.complete == False:
            print(i)
            self.update_Jacobi()
            # self.update_Gauss_Seidel()
            self.check_accuracy()
            i +=1
        print("Accuracy of %.6f reached"  %(self.tolerance))

        np.save("phi_wire",self.phi)
        return i

    def load(self):
        self.phi = np.load("phi.npy")
#------------------------------------------------------------------------------
    def fit_phi(self):
        phi_x = self.phi[int(self.N/2):,int(self.N/2),int(self.N/2)]
        x = range(0,int(self.N/2))

        x = np.log(x)
        y = np.log(phi_x)

        # coefficients = np.polyfit(y, x, 1)
        # polynomial = np.poly1d(coefficients)
        # ys = polynomial(b)
        # plt.plot(y, x)
        # plt.plot(b, ys)

        # np.polyfit(x,y,1)
        plt.xlabel("Log of distance from charge")
        plt.ylabel("Log of Electrostatic Potential")
        plt.title("Log log plot of Phi as function of distance from charge")
        plt.scatter(x,y)
        plt.show()

        # return np.polyfit(x,y,4)

    def fit_E(self):
        E_x = abs(self.Ex[int(self.N/2) + 1:,int(self.N/2),int(self.N/2)])
        print (E_x)
        x = range(1,int(self.N/2))

        x = np.log(x)
        y = np.log(E_x)

        np.polyfit(x,y,1)
        plt.xlabel("Log of distance from charge")
        plt.ylabel("Log of Electric Field Strength")
        plt.title("Log log plot of E as function of distance from charge")
        plt.scatter(x,y)
        plt.show()

        # return np.polyfit(x,y,4)

    def fit_B(self):
        B_x = abs(self.Bx[int(self.N/2) + 1:,int(self.N/2)])
        print (B_x)
        x = range(1,int(self.N/2))

        x = np.log(x)
        y = np.log(B_x)

        np.polyfit(x,y,1)
        plt.xlabel("Log of distance from charge")
        plt.ylabel("Log of Magnetic Field Strength")
        plt.title("Log log plot of B as function of distance from charge")
        plt.scatter(x,y)
        plt.show()

        return np.polyfit(x[:5],y[:5],1)

#------------------------------------------------------------------------------
def omega_test():
    convergence_times = []
    omega = np.linspace(1,2,10)[0:-1]
    for i in omega:
        x = Poisson(50)
        x.place_central_charge()
        print (i)
        x.omega = i
        t = x.run()
        convergence_times.append(t)
        x.complete = False
    np.savetxt("convergence_times_new.txt", convergence_times)
    np.savetxt("omega_vals_new.txt", omega)
    plt.plot(omega,convergence_times)
    plt.show()

def main():
    N = sys.argv[1]
    Accuracy = sys.argv[2]
    x = Poisson(N,Accuracy)
    # x.place_central_charge()
    x.place_wire()

    x.run()

    # x.load()

    x.get_E()
    x.get_B()


    # x.plot_electrostatic_potential_x()
    # x.plot_electric_field_x()

    # x.B_field_plot()

    # x.norm_fields()

    # x.E_vector_plot_XY()

    # x.B_vector_plot_XY()


    # x.fit_phi()
    # x.fit_E()
    print (x.fit_B())


main()
# omega_test()
