import pymbar
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

interface_list = [-0.625, -0.600, -0.580, -0.560, -0.540, -0.520, -0.500, -0.425, -0.350]
interface_list = np.array(interface_list)
K_int = len(interface_list)

theta_ik = []
 
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('xtick', labelsize=28)
plt.rc('ytick', labelsize=28)

for k, biasval in enumerate(interface_list):
    data = np.genfromtxt('theta{:1.3f}.txt'.format(biasval))
#     data = data.reshape((-1, 240))
#     data_t = np.mean(data, axis=1)
    theta_ik.append(data[:])
theta_ik = np.array(theta_ik)

color_list=iter(plt.cm.rainbow(np.linspace(0,1,K_int)))

bin_centers = np.linspace(-1.0, 0.5, 250)

for index in range(K_int):
    c = next(color_list)
    plt.hist(theta_ik[index, :], bins=bin_centers, normed=True, label=r'$\theta = {}$'.format(interface_list[index]), histtype='stepfilled', alpha=0.5, color=c)
    plt.xlabel(r'$\theta$', fontsize=28)
    plt.ylabel(r'$P(\theta)$', fontsize=28)
    plt.legend(loc='best', fontsize=20)
    plt.show()



