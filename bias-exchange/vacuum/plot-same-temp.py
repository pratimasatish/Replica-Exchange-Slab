import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

lwidth = 4.0
plt.rc('text', usetex=True, fontsize=30)
rcParams['text.latex.preamble'] = [r'\usepackage{helvet} \usepackage{sfmath}', r'\usepackage{upgreek}' ]
rcParams['axes.linewidth'] = 1.0*lwidth
rcParams['xtick.major.width'] = 1.0*lwidth
rcParams['xtick.major.size']  = 2.0*lwidth
rcParams['ytick.major.width'] = 1.0*lwidth
rcParams['ytick.major.size']  = 2.0*lwidth
plt.rc('lines', linewidth=4)
plt.rc('legend', frameon=False)

data_370 = np.genfromtxt('f_i-370.txt')
data_380 = np.genfromtxt('f_i-380.txt')

plt.figure(figsize=(16,9))
plt.plot(data_370[:,0], data_370[:,1], '#7fbf7b', marker='^', markersize=10, label='370K')
plt.plot(data_370[:,0], data_370[:,1], '#7fbf7b', lw=4, alpha=0.7)
plt.plot(data_380[1:,0], data_380[1:,1], color='#af8dc3', marker='s', markersize=10, label='380K')
plt.plot(data_380[1:,0], data_380[1:,1], '#af8dc3', lw=4, alpha=0.7)
leg = plt.legend(loc='upper left')
for legobj in leg.legendHandles:
    legobj.set_linewidth(3.5)
plt.xlabel(r'$\Theta_z\,\mathsf{(rad)}$')
plt.ylabel(r'$\mathsf{free\,\,energy\,\,density\,\,(kcal/mol}\mbox{-}\mathsf{rad)}$')
plt.savefig('free-en-370+380K.pdf')

# plt.figure()
# plt.plot(data_370[:,0], data_370[:,2], 'r^', markersize=8, label='370K')
# plt.plot(data_380[:,0], data_380[:,2], 'gs', markersize=8, label='380K')
# plt.legend(loc='upper left', fontsize=28)
# plt.xlabel(r'\boldsymbol{$\langle\theta_z\rangle}', fontsize=28)
# plt.ylabel(r'$\textbf{energy density (kcal/mol-rad)}$', fontsize=28)
# plt.show()

plt.figure(figsize=(12,9))
plt.plot(data_370[:,0], data_370[:,3], color='#7fbf7b', marker='^', markersize=10, label='370K')
plt.plot(data_380[:,0], data_380[:,3], color='#af8dc3', marker='s', markersize=10, label='380K')
plt.legend(loc='upper left')
plt.xlabel(r'$\Theta_z\,\mathsf{(rad)}$')
plt.ylabel(r'$\mathsf{entropy\,\,density\,\,(kcal/mol}\mbox{-}\mathsf{K}\mbox{-}\mathsf{rad)}$')
plt.xticks(np.arange(-0.8, 0, 0.1))
plt.yticks(np.arange(0, 4.5, 0.5))
leg = plt.legend(loc='upper left')
for legobj in leg.legendHandles:
    legobj.set_linewidth(3.5)
plt.savefig('entropies-370+380K.pdf')

plt.figure(figsize=(14,9))
plt.plot(data_370[:,0], data_370[:,2], color='#7fbf7b', marker='^', markersize=10, label='370K')
plt.plot(data_380[:,0], data_380[:,2], color='#af8dc3', marker='s', markersize=10, label='380K')
plt.legend(loc='upper left')
plt.xlabel(r'$\Theta_z\,\mathsf{(rad)}$')
plt.ylabel(r'$\mathsf{energy\,\,density\,\,(kcal/mol}\mbox{-}\mathsf{K}\mbox{-}\mathsf{rad)}$')
plt.xticks(np.arange(-0.8, 0, 0.1))
# plt.yticks(np.arange(0, 4.5, 0.5))
leg = plt.legend(loc='upper left')
for legobj in leg.legendHandles:
    legobj.set_linewidth(3.5)
plt.savefig('energies-370+380K.pdf')


data = np.genfromtxt('f_i-combined-375.6.txt')
temp = 375.6

plt.figure(figsize=(16,9))
# plt.plot(data[:,0], data[:,1], 'k^', markersize=10, label='{}K'.format(temp))
plt.plot(data[:,0], data[:,1], 'k', lw=6, label='{}K'.format(temp))
# plt.legend(loc='best', fontsize=28)
plt.xlabel(r'$\Theta_z\,\mathsf{(rad)}$')
plt.ylabel(r'$\mathsf{free\,\,energy\,\,density\,\,(kcal/mol}\mbox{-}\mathsf{rad)}$')
plt.xlim(-0.8, 0.0)
#plt.tight_layout()
plt.savefig('free-en-coex-375.6K.pdf')
