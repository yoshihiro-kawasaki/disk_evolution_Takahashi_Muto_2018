import numpy as np 
import matplotlib.pyplot as plt 

file = "../be_sphere.txt"
data = np.loadtxt(file, dtype=float, skiprows=1).T

x = data[0]       # normalized radius
psi = data[3]     # normalized gravitatio potential
dpsi = data[4]    
del data 
rho = np.exp(-psi) # normalized mass density

plt.figure(figsize=(12, 10))
plt.plot(x, psi, label="psi")
plt.plot(x, rho, label="rho")
plt.tick_params(labelsize=30)
plt.legend(fontsize=30)
plt.show()