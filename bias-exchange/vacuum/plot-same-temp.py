import numpy as np
import matplotlib.pyplot as plt

# data_370 = np.genfromtxt('f_i-370.txt')
# data_380 = np.genfromtxt('f_i-380.txt')
# 
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')
# plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
# plt.rc('font', size=28)
# plt.rc('xtick', labelsize=28)
# plt.rc('ytick', labelsize=28)
# 
# plt.figure()
# plt.plot(data_370[:,0], data_370[:,1], 'r^', markersize=8, label='370K')
# plt.plot(data_380[:,0], data_380[:,1], 'gs', markersize=8, label='380K')
# plt.legend(loc='upper left', fontsize=28)
# plt.xlabel(r'\boldsymbol{$\langle\theta_z\rangle}', fontsize=28)
# plt.ylabel(r'$\textbf{free energy density (kcal/mol-rad)}$', fontsize=28)
# plt.show()
# 
# plt.figure()
# plt.plot(data_370[:,0], data_370[:,2], 'r^', markersize=8, label='370K')
# plt.plot(data_380[:,0], data_380[:,2], 'gs', markersize=8, label='380K')
# plt.legend(loc='upper left', fontsize=28)
# plt.xlabel(r'\boldsymbol{$\langle\theta_z\rangle}', fontsize=28)
# plt.ylabel(r'$\textbf{energy density (kcal/mol-rad)}$', fontsize=28)
# plt.show()
# 
# plt.figure()
# plt.plot(data_370[:,0], data_370[:,3], 'r^', markersize=8, label='370K')
# plt.plot(data_380[:,0], data_380[:,3], 'gs', markersize=8, label='380K')
# plt.legend(loc='upper left', fontsize=28)
# plt.xlabel(r'\boldsymbol{$\langle\theta_z\rangle}', fontsize=28)
# plt.ylabel(r'$\textbf{entropy density (kcal/mol-K-rad)}$', fontsize=28)
# plt.show()


data = np.genfromtxt('f_i-combined-375.6.txt')
temp = 375.6

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
plt.rc('font', size=28)
plt.rc('xtick', labelsize=28)
plt.rc('ytick', labelsize=28)

plt.figure()
plt.plot(data[:,0], data[:,1], 'r^', markersize=8, label='{}K'.format(temp))
plt.legend(loc='upper left', fontsize=28)
plt.xlabel(r'\boldsymbol{$\langle\theta_z\rangle}', fontsize=28)
plt.ylabel(r'$\textbf{free energy density (kcal/mol-rad)}$', fontsize=28)
plt.xlim(-0.8, 0.0)
plt.show()
