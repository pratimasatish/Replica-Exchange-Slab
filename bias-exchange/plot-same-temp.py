import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

data_365 = np.genfromtxt('f_i-365.txt')
data_375 = np.genfromtxt('f_i-375.txt')

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

plt.figure(figsize=(16,9))
plt.plot(data_365[:,0], data_365[:,1], '#7fbf7b', marker='^', markersize=10, label='365K')
plt.plot(data_365[:,0], data_365[:,1], '#7fbf7b', lw=4, alpha=0.7)
plt.plot(data_375[:,0], data_375[:,1], '#af8dc3', marker='s', markersize=8, label='375K')
plt.plot(data_375[:,0], data_375[:,1], '#af8dc3', lw=4, alpha=0.7)
leg = plt.legend(loc='upper right')
for legobj in leg.legendHandles:
    legobj.set_linewidth(3.5)
plt.xlabel(r'$\Theta_z\,\mathsf{(rad)}$')
plt.ylabel(r'$\mathsf{free\,\,energy\,\,density\,\,(kcal/mol}\mbox{-}\mathsf{rad)}$')
plt.savefig('free-en-365+375K.pdf')

plt.figure(figsize=(12,9))
plt.plot(data_365[:,0], data_365[:,3], color='#7fbf7b', marker='^', markersize=10, label='365K')
plt.plot(data_375[:,0], data_375[:,3], color='#af8dc3', marker='s', markersize=10, label='375K')
plt.legend(loc='upper left')
plt.xlabel(r'$\langle\theta_z\rangle$')
plt.ylabel(r'$\mathsf{entropy\,\,density\,\,(kcal/mol}\mbox{-}\mathsf{K}\mbox{-}\mathsf{rad)}$')
plt.xticks(np.arange(-0.8, 0, 0.1))
plt.yticks(np.arange(0, 4.5, 0.5))
leg = plt.legend(loc='upper left')
for legobj in leg.legendHandles:
    legobj.set_linewidth(3.5)
plt.savefig('entropies-365+375K.pdf')

plt.figure(figsize=(14,9))
plt.plot(data_365[:,0], data_365[:,2], color='#7fbf7b', marker='^', markersize=10, label='365K')
plt.plot(data_375[:,0], data_375[:,2], color='#af8dc3', marker='s', markersize=10, label='375K')
plt.legend(loc='upper left')
plt.xlabel(r'$\Theta_z\,\mathsf{(rad)}$')
plt.ylabel(r'$\mathsf{energy\,\,density\,\,(kcal/mol}\mbox{-}\mathsf{rad)}$')
plt.xticks(np.arange(-0.8, 0, 0.1))
# plt.yticks(np.arange(0, 4.5, 0.5))
leg = plt.legend(loc='upper left')
for legobj in leg.legendHandles:
    legobj.set_linewidth(3.5)
plt.savefig('energies-365+375K.pdf')


# plt.figure()
# plt.plot(data_365[:,0], data_365[:,2], 'r^', markersize=8, label='365K')
# plt.plot(data_375[:,0], data_375[:,2], 'gs', markersize=8, label='375K')
# plt.legend(loc='upper left', fontsize=28)
# plt.xlabel(r'\boldsymbol{$\langle\theta_z\rangle}', fontsize=28)
# plt.ylabel(r'$\textbf{energy density (kcal/mol-rad)}$', fontsize=28)
# plt.show()

data = np.genfromtxt('f_i-368-from-365.txt')
temp = 368

plt.figure(figsize=(16,9))
# plt.plot(data[:,0], data[:,1], 'k^', markersize=10, label='{}K'.format(temp))
plt.plot(data[:,0], data[:,1], 'k', lw=6, label='{}K'.format(temp))
# plt.legend(loc='best', fontsize=28)
plt.xlabel(r'$\Theta_z\,\mathsf{(rad)}$')
plt.ylabel(r'$\mathsf{free\,\,energy\,\,density\,\,(kcal/mol}\mbox{-}\mathsf{rad)}$')
plt.xlim(-0.8, 0.0)
#plt.tight_layout()
plt.savefig('free-en-coex-368K.pdf')


# data = np.genfromtxt('f_i-368-from-365.txt')
# temp = 368.1
# 
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')
# plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
# plt.rc('font', size=28)
# plt.rc('xtick', labelsize=28)
# plt.rc('ytick', labelsize=28)
# 
# plt.figure()
# plt.plot(data[:,0], data[:,1], 'r^', markersize=8, label='{}K'.format(temp))
# plt.legend(loc='upper left', fontsize=28)
# plt.xlabel(r'\boldsymbol{$\langle\theta_z\rangle}', fontsize=28)
# plt.ylabel(r'$\textbf{free energy density (kcal/mol-rad)}$', fontsize=28)
# plt.xlim(-0.8, 0.0)
# plt.show()
