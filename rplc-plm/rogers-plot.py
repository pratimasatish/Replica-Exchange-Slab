import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('log.stripped')

plt.figure()
plt.plot(data[:,0], data[:,1], 'r')
plt.plot(data[:,0], data[:,2], 'b')
plt.ylim(-0.5, 1.5)
plt.show()

