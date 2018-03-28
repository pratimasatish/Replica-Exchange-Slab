import numpy as np
import matplotlib.pyplot as plt

data_350 = np.genfromtxt('350K/364-from-350.txt')
data_365 = np.genfromtxt('365K/370-from-365.txt')
data_370 = np.genfromtxt('370K/372-from-370.txt')
data_375 = np.genfromtxt('375K/373-from-375.txt')

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('xtick', labelsize=28)
plt.rc('ytick', labelsize=28)

plt.figure()
plt.plot(data_350[:,0], data_350[:,1], 'rv', markersize=8, label='from 350K: transition at 364K')
plt.plot(data_365[:,0], data_365[:,1], 'bo', markersize=8, label='from 365K: transition at 370K')
plt.plot(data_370[:,0], data_370[:,1], 'k^', markersize=8, label='from 370K: transition at 372K')
plt.plot(data_375[:,0], data_375[:,1], 'gs', markersize=8, label='from 375K: transition at 373K')
plt.legend(loc='upper left')
plt.xlabel(r'$\langle\theta_z\rangle', fontsize=28)
plt.ylabel(r'$-\beta F', fontsize=28)
plt.show()



