import numpy as np
import matplotlib.pyplot as plt

x = np.arange(1,3,0.1)

Glauber_Mean_Energy = np.loadtxt("Data/Mean_Energy_Glauber.txt")
Glauber_Mean_Magnetiation = np.loadtxt("Data/Mean_Magnetisation_Glauber.txt")
Kawasaki_Mean_Energy = np.loadtxt("Data/Mean_Energy_Kawasaki.txt")
Glauber_Heat_Capacity = np.loadtxt("Data/GlauberHeatCapacity.txt")
Glauber_Susceptibility = np.loadtxt("Data/GlauberSuscepibility.txt")
Kawasaki_Heat_Capacity = np.loadtxt("Data/KawasakiHeatCapacity.txt")

Kawasaki_C_Error = np.loadtxt("Data/KawasakiC_Error.txt")
Glauber_C_Error = np.loadtxt("Data/GlauberC_Error.txt")
Glauber_Chi_Error = np.loadtxt("Data/GlauberChi_Error.txt")


plt.errorbar(x, Kawasaki_Heat_Capacity, yerr=Kawasaki_C_Error, color = "Black")
plt.title("Kawasaki Heat Capacity")
plt.ylabel("Heat Capacity")
plt.xlabel("Temperature (K)")
x1,x2,y1,y2 = plt.axis()
plt.axis((1,x2,y1,y2))
plt.savefig("Kawasaki Heat Capacity.png")
plt.show()

plt.errorbar(x, Glauber_Heat_Capacity, yerr=Glauber_C_Error, color = "Black")
plt.title("Glauber Heat Capacity")
plt.ylabel("Heat Capacity")
plt.xlabel("Temperature (K)")
x1,x2,y1,y2 = plt.axis()
plt.axis((1,x2,y1,y2))
plt.savefig("Glauber Heat Capacity.png")

plt.show()

plt.errorbar(x, Glauber_Susceptibility, yerr=Glauber_Chi_Error, color = "Black")
plt.title("Glauber Suscepibility")
plt.ylabel("Suscepibility")
plt.xlabel("Temperature (K)")
x1,x2,y1,y2 = plt.axis()
plt.axis((1,x2,y1,y2))
plt.savefig("Glauber Suscepibility.png")
plt.show()

plt.plot(x,Glauber_Mean_Energy, color = "Black")
plt.title("Glauber Mean Energy")
plt.ylabel("Energy")
plt.xlabel("Temperature (K)")
x1,x2,y1,y2 = plt.axis()
plt.axis((1,x2,y1,y2))
plt.savefig("Glauber Mean Energy.png")
plt.show()

plt.plot(x,Glauber_Mean_Magnetiation, color = "Black")
plt.title("Glauber Mean Magnetisation")
plt.ylabel("Magnetisation")
plt.xlabel("Temperature (K)")
x1,x2,y1,y2 = plt.axis()
plt.axis((1,x2,y1,y2))
plt.savefig("Glauber Mean Magnetisation.png")
plt.show()

plt.plot(x,Kawasaki_Mean_Energy, color = "Black")
plt.title("Kawasaki Mean Energy")
plt.ylabel("Energy")
plt.xlabel("Temperature (K)")
x1,x2,y1,y2 = plt.axis()
plt.axis((1,x2,y1,y2))
plt.savefig("Kawasaki Mean Energy")
plt.show()