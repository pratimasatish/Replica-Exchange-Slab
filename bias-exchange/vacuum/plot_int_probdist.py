import pymbar
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

namelist = ['ord', 'disord']
K_int = len(namelist)

theta_ik = []
 
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('xtick', labelsize=28)
plt.rc('ytick', labelsize=28)

for biasval in namelist:
    data = np.genfromtxt('theta' + biasval + '.txt')
#     data = data.reshape((-1, 240))
#     data_t = np.mean(data, axis=1)
    theta_ik.append(data[:])
theta_ik = np.array(theta_ik)

color_list=iter(plt.cm.rainbow(np.linspace(0,1,K_int)))

bin_centers = np.linspace(-1.0, 1.0, 250)

for index in range(K_int):
#     c = next(color_list)
    if index == 0:
        plt.hist(theta_ik[index, :], bins=bin_centers, normed=True, histtype='stepfilled', alpha=0.7, color='r')
    else:
        plt.hist(theta_ik[index, :], bins=bin_centers, normed=True, histtype='stepfilled', alpha=0.7, color='b')
    plt.xlabel(r'$\theta$', fontsize=28)
    plt.ylabel(r'$P(\theta)$', fontsize=28)
#     plt.legend(loc='best', fontsize=20)
    plt.show()



