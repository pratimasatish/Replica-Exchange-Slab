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

data = np.genfromtxt('etadiff-vs-temp.txt')

new_data = []
for i in range(len(data)):
    if (data[i,3] == 1):
        new_data.append(data[i,1])
    elif (data[i,3] == 0):
        new_data.append(data[i,1]*0.5 + data[i,2]*0.5)
    else:
        new_data.append(data[i,2])
new_data = np.array(new_data) 

plt.figure(figsize=(14,10))
plt.plot(data[:,0], data[:,1], color='#7fbf7b', marker='^', markersize=12, lw=3, label=r'$\langle\theta_\mathsf{z}\rangle_\mathsf{ord}$')
plt.plot(data[:,0], data[:,2], color='#af8dc3', marker='s', markersize=12, lw=3, label=r'$\langle\theta_\mathsf{z}\rangle_\mathsf{disord}$')
plt.plot(data[:,0], data[:,2]-data[:,1], color='#253494', marker='D', markersize=12, lw=4, alpha=0.8, label=r'$\Delta\langle\theta_{\mathsf{z}}\rangle$')
# plt.plot(data[:,0], data[:,2]-data[:,1], color='#253494', linewidth=4, alpha=0.7)
plt.hlines(np.mean(data[:,2]-data[:,1]), 363, 377, 'k', linewidth=3, linestyle='dashed', label=r'$\Delta\langle\theta_\mathsf{z}\rangle_\mathsf{avg}$')
plt.legend(loc=(0.6, 0.5))
# plt.legend(loc='best')
plt.xlabel(r'$\mathsf{T (K)}$')
plt.ylabel(r'$\langle\theta_{\mathsf{z}}\rangle$')
plt.xlim(363, 377)
plt.ylim(-0.8, 0.6)
plt.savefig('etadiff-vs-T-solv.pdf')

# plt.figure(figsize=(12,9))
# plt.plot(data[:,0], new_data, color='#253494', marker='^', markersize=12)
# plt.plot(data[:,0], new_data, color='#253494', linewidth=4, alpha=0.7)
# plt.hlines(np.mean(data[:,2]-data[:,1]), 368, 382, 'k', linewidth=3, linestyle='dashed')
# # plt.legend(loc='upper left', fontsize=28)
# plt.xlabel(r'$\textrm{T (K)}$', fontsize=32)
# plt.ylabel(r'$\langle\theta_{\textrm{z}}\rangle$', fontsize=32)
# plt.xlim(368, 382)
# plt.ylim(-0.8, -0.2)
# plt.savefig('etadiff-vs-T.pdf')
# 



