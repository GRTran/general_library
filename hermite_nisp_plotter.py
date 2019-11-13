import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('random_inputs.dat')

plt.scatter(data[:,0], data[:,1])
plt.show()

plt.scatter(data[:,0],data[:,2])
plt.scatter(data[:,0],data[:,3])
plt.show()
