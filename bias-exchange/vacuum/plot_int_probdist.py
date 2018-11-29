import pymbar
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

namelist = ['ord', 'disord']
K_int = len(namelist)

theta_ik = []
 
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
plt.rc('xtick', labelsize=28)
plt.rc('ytick', labelsize=28)

for biasval in namelist:
    data = np.genfromtxt('theta' + biasval + '.txt')
#     data = data.reshape((-1, 240))
#     data_t = np.mean(data, axis=1)
    theta_ik.append(data[:])
theta_ik = np.array(theta_ik)
bin_centers = np.linspace(-1.0, 1.0, 250)

# color_list=iter(plt.cm.rainbow(np.linspace(0,1,K_int)))
# 
# for index in range(K_int):
# #     c = next(color_list)
#     if index == 0:
#         plt.hist(theta_ik[index, :], bins=bin_centers, normed=True, histtype='stepfilled', alpha=0.7, color='r')
#     else:
#         plt.hist(theta_ik[index, :], bins=bin_centers, normed=True, histtype='stepfilled', alpha=0.7, color='b')
#     plt.xlabel(r'$\theta$', fontsize=28)
#     plt.ylabel(r'$P(\theta)$', fontsize=28)
# #     plt.legend(loc='best', fontsize=20)
#     plt.show()

theta_ord_red = theta_ik[0][np.where(theta_ik[0]<-0.77778)]
theta_ord_blue = theta_ik[0][np.where(theta_ik[0]>=-0.77778)]
theta_disord_red = theta_ik[1][np.where(theta_ik[1]<-0.77778)]
theta_disord_blue = theta_ik[1][np.where(theta_ik[1]>=-0.77778)]

plt.hist(theta_ord_red, bins=bin_centers, histtype='stepfilled', alpha=0.7, color='r')
plt.hist(theta_disord_red, bins=bin_centers, histtype='stepfilled', alpha=0.7, color='r')
plt.hist(theta_ord_blue, bins=bin_centers, histtype='stepfilled', alpha=0.7, color='b')
plt.hist(theta_disord_blue, bins=bin_centers, histtype='stepfilled', alpha=0.7, color='b')
plt.xlabel(r'$\boldsymbol{\theta_z}$', fontsize=28)
plt.ylabel(r'$\boldsymbol{P(\theta_z)}$', fontsize=28)
plt.show()

theta_ord_green = theta_ik[0][np.where(theta_ik[0]<-0.575)]
theta_ord_purple = theta_ik[0][np.where(theta_ik[0]>=-0.575)]
theta_disord_green = theta_ik[1][np.where(theta_ik[1]<-0.575)]
theta_disord_purple = theta_ik[1][np.where(theta_ik[1]>=-0.575)]

plt.hist(theta_ord_green, bins=bin_centers, histtype='stepfilled', alpha=0.7, color='g')
plt.hist(theta_disord_green, bins=bin_centers, histtype='stepfilled', alpha=0.7, color='g')
plt.hist(theta_ord_purple, bins=bin_centers, histtype='stepfilled', alpha=0.7, color='purple')
plt.hist(theta_disord_purple, bins=bin_centers, histtype='stepfilled', alpha=0.7, color='purple')
plt.xlabel(r'$\boldsymbol{\theta_z}$', fontsize=28)
plt.ylabel(r'$\boldsymbol{P(\theta_z)}$', fontsize=28)
plt.show()


