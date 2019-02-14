import numpy as np

Kawasaki_Energy = np.loadtxt("Data/KawasakiEnergy.txt")
Glauber_Magnetisation = np.loadtxt("Data/GlauberMagnetisation.txt")
Glauber_Energy = np.loadtxt("Data/GlauberEnergy.txt")

def get_mean_vals(Data):
  mean_vals = []
  for i in range(len(Data)):
    mean_vals.append(np.mean(Data[i]))

  return mean_vals

np.savetxt(str("Data/Mean_Energy_Kawasaki.txt"), get_mean_vals(Kawasaki_Energy))
np.savetxt(str("Data/Mean_Energy_Glauber.txt"), get_mean_vals(Glauber_Energy))
np.savetxt(str("Data/Mean_Magnetisation_Glauber.txt"), get_mean_vals(Glauber_Magnetisation))
