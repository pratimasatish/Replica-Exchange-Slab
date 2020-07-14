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

bin_edges = np.linspace(-1.0, 1.0, 201)
bins = 0.5 * (bin_edges[1:] + bin_edges[:-1])

data_DO = np.genfromtxt('theta375.txt')
data_O = np.genfromtxt('theta360.txt')

bin_edges = np.linspace(-1.0, 1.0, 201)
bins = 0.5 * (bin_edges[1:] + bin_edges[:-1])

hist_DO,_ = np.histogram(data_DO, bins=bin_edges, density=True)
hist_O,_ = np.histogram(data_O, bins=bin_edges, density=True)

plt.figure(figsize=(12,9))
plt.plot(bins, hist_O, '#7fbf7b', lw=6, label='Ordered phase')
plt.plot(bins, hist_DO, '#af8dc3', lw=6, label='Disordered phase')
plt.xlabel(r'$\theta_z$')
plt.ylabel(r'$P(\theta_z)$')
plt.legend(loc='best')
plt.ylim(0, 6)
# plt.show()
plt.savefig('compare-solv-pdist-thz.pdf')

tz_DO = np.reshape(data_DO, (-1, 240))
data_DO = np.mean(tz_DO, axis=1)
tz_O = np.reshape(data_O, (-1, 240))
data_O = np.mean(tz_O, axis=1)

bin_DO_edges = np.linspace(min(data_DO), max(data_DO), 31)
bin_DO = 0.5 * (bin_DO_edges[1:] + bin_DO_edges[:-1]) 
print(min(data_O), max(data_O))
bin_O_edges = np.linspace(min(data_O), max(data_O), 31)
bin_O = 0.5 * (bin_O_edges[1:] + bin_O_edges[:-1]) 

hist_DO,_ = np.histogram(data_DO, bins=bin_DO_edges, density=True)
hist_O,_ = np.histogram(data_O, bins=bin_O_edges, density=True)

plt.figure(figsize=(12,9))
plt.plot(bin_O, hist_O, '#7fbf7b', lw=6, label='Ordered phase')
plt.plot(bin_DO, hist_DO, '#af8dc3', lw=6, label='Disordered phase')
plt.xlabel(r'$\Theta_z$')
plt.ylabel(r'$P(\Theta_z)$')
plt.legend(loc='best')
plt.ylim(0, 65)
plt.xlim(-0.8, 0)
# plt.show()
plt.savefig('compare-solv-pdist-avethz.pdf')

