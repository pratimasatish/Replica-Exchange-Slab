import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('log.stripped')

plt.figure()
plt.plot(data[:,0], data[:,1], 'rx')
plt.plot(data[:,0], data[:,2], 'b')
plt.plot(data[:,0], data[:,3], 'k')
plt.plot(data[:,0], data[:,4], 'g')
plt.plot(data[:,0], data[:,5], 'y')
plt.plot(data[:,0], data[:,6], 'o')
plt.ylim(-0.5, 6.5)
plt.show()

