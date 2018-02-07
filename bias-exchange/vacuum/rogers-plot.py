import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('log.stripped')

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('xtick', labelsize=28)
plt.rc('ytick', labelsize=28)

plt.figure()
plt.plot(data[:,0], data[:,1], 'rx')
plt.plot(data[:,0], data[:,2], 'b')
plt.plot(data[:,0], data[:,3], 'k')
plt.plot(data[:,0], data[:,4], 'g')
plt.plot(data[:,0], data[:,5], 'y')
plt.plot(data[:,0], data[:,6], 'o')
plt.plot(data[:,0], data[:,7], 'm')
plt.plot(data[:,0], data[:,8], 'c')
plt.plot(data[:,0], data[:,9], '#eeefff')
plt.plot(data[:,0], data[:,10], '0.75')
plt.plot(data[:,0], data[:,11], '0.85')
plt.plot(data[:,0], data[:,12], '0.35')
plt.plot(data[:,0], data[:,13], '0.25')
plt.plot(data[:,0], data[:,14], '0.65')
plt.plot(data[:,0], data[:,15], '0.95')
plt.ylim(-0.5, 15.5)
plt.xlabel("time", fontsize=28)
plt.ylabel("replica id", fontsize=28)
plt.show()

