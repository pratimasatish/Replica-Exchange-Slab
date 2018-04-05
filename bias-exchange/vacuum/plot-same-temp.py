import numpy as np
import matplotlib.pyplot as plt

data_370 = np.genfromtxt('370K/372.5-from-370.txt')
data_375 = np.genfromtxt('375K/372.5-from-375.txt')
data_trans = np.genfromtxt('375K/373-from-375.txt')

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('xtick', labelsize=28)
plt.rc('ytick', labelsize=28)

plt.figure()
plt.plot(data_370[:,0], data_370[:,1], 'k^', markersize=8, label='from 370K: at 372.5K')
plt.plot(data_375[:,0], data_375[:,1], 'gs', markersize=8, label='from 375K: at 372.5K')
plt.plot(data_trans[:,0], data_trans[:,1], 'ro', markersize=8, label='predicted coexistence free energy')
plt.legend(loc='upper left')
plt.xlabel(r'$\langle\theta_z\rangle', fontsize=28)
plt.ylabel(r'$-\beta F', fontsize=28)
plt.show()



