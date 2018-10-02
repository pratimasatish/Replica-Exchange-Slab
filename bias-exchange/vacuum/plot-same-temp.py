import numpy as np
import matplotlib.pyplot as plt

data_370 = np.genfromtxt('f_i-370.txt')
data_380 = np.genfromtxt('f_i-380.txt')

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('xtick', labelsize=28)
plt.rc('ytick', labelsize=28)

plt.figure()
plt.plot(data_370[:,0], data_370[:,1], 'r^', markersize=8, label='370K')
plt.plot(data_380[:,0], data_380[:,1], 'gs', markersize=8, label='380K')
plt.legend(loc='upper left')
plt.xlabel(r'$\langle\theta_z\rangle', fontsize=28)
plt.ylabel(r'$F \textrm{ (kcal/mol)}$', fontsize=28)
plt.show()

plt.figure()
plt.plot(data_370[:,0], data_370[:,2], 'r^', markersize=8, label='370K')
plt.plot(data_380[:,0], data_380[:,2], 'gs', markersize=8, label='380K')
plt.legend(loc='upper left')
plt.xlabel(r'$\langle\theta_z\rangle', fontsize=28)
plt.ylabel(r'$E \textrm{ (kcal/mol)}$', fontsize=28)
plt.show()

plt.figure()
plt.plot(data_370[:,0], data_370[:,3], 'r^', markersize=8, label='370K')
plt.plot(data_380[:,0], data_380[:,3], 'gs', markersize=8, label='380K')
plt.legend(loc='upper left')
plt.xlabel(r'$\langle\theta_z\rangle', fontsize=28)
plt.ylabel(r'$S \textrm{ (kcal/mol-K)}$', fontsize=28)
plt.show()

