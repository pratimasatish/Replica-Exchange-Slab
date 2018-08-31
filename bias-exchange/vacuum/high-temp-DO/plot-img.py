import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('tcf-noavg-removed-1000K.txt')

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=28)
plt.rc('font', weight='bold')
plt.rc('xtick', labelsize=28)
plt.rc('ytick', labelsize=28)

plt.plot(data[:,0], data[:,1], 'ro')
plt.plot(data[:,0], data[:,1], 'r', linewidth=4, alpha=0.5)
plt.show()

