import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('etadiff-vs-temp.txt')

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
plt.rc('font', size=28)
plt.rc('xtick', labelsize=28)
plt.rc('ytick', labelsize=28)

plt.figure()
plt.plot(data[:,0], (data[:,2]-data[:,1]), 'r^', markersize=12)
plt.plot(data[:,0], (data[:,2]-data[:,1]), 'r', linewidth=4, alpha=0.7)
plt.hlines(np.mean(data[:,2]-data[:,1]), 368, 382, 'k', linewidth=3, linestyle='dashed')
# plt.legend(loc='upper left', fontsize=28)
plt.xlabel(r'$\textbf{T (K)}$', fontsize=32)
plt.ylabel(r'$\boldsymbol{\Delta\theta_z}$', fontsize=32)
plt.xlim(368, 382)
plt.ylim(0.400, 0.430)
plt.show()
