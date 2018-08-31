import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('xtick', labelsize=28)
plt.rc('ytick', labelsize=28)

lotemp = np.genfromtxt('g_r_365.dat')
hitemp = np.genfromtxt('g_r_1000.dat')

kB = 1.3806503 * 6.0221415 / 4184.0
lobeta = 1/(kB * 365)
hibeta = 1/(kB * 1000)

# lotemp_tot = np.sum(lotemp[:,1])
# hitemp_tot = np.sum(hitemp[:,1])

lotemp_tot = np.mean(lotemp[:,1][np.where(lotemp[:,0] >= 6)])
hitemp_tot = np.mean(hitemp[:,1][np.where(hitemp[:,0] >= 6)])

lotemp[:,1] /= lotemp_tot
hitemp[:,1] /= hitemp_tot

eps = 0.0914
sig = 3.95
lj_pot = np.zeros(len(lotemp[:,0]))
for i in range(len(lotemp[:,0])):
    lj_pot[i] = 4 * eps * ((sig / lotemp[:,0][i])**12 - (sig / lotemp[:,0][i])**6)

plt.plot(lotemp[:,0], -np.log(lotemp[:,1]) / lobeta, 'ro', label='365K')
plt.plot(lotemp[:,0], -np.log(lotemp[:,1]) / lobeta, 'r', linewidth=4, alpha=0.4)
plt.plot(hitemp[:,0], -np.log(hitemp[:,1]) / hibeta, 'bv', label='1000K')
plt.plot(hitemp[:,0], -np.log(hitemp[:,1]) / hibeta, 'b', linewidth=4, alpha=0.4)
plt.plot(lotemp[:,0], lj_pot, 'ks', label='LJ pot')
plt.plot(lotemp[:,0], lj_pot, 'k', linewidth=4, alpha=0.4)
plt.legend(loc='upper right')
plt.show()

plt.plot(lotemp[:,0], lotemp[:,1], 'ro', label='365K')
plt.plot(lotemp[:,0], lotemp[:,1], 'r', linewidth=4, alpha=0.4)
plt.plot(hitemp[:,0], hitemp[:,1], 'bv', label='1000K')
plt.plot(hitemp[:,0], hitemp[:,1], 'b', linewidth=4, alpha=0.4)
plt.legend(loc='upper right')
plt.show()

