import numpy as np
import argparse
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math

parser = argparse.ArgumentParser(description="")
parser.add_argument("-bias", type=str, help="bias value to analyse")
parser.add_argument("-steps", type=int, default=5000, help="number of time steps to analyse")
parser.add_argument("-pbc", action='store_true', help="whether to apply PBCs or not")
parser.add_argument("-savefigs", action='store_true', help="whether to save figures of CG lattice to make movie or not")
# parser.add_argument("-remove_NN", action='store_true', help="whether to remove sites with less than 4 NN's")
args = parser.parse_args()

data = np.genfromtxt('bottheta' + args.bias + '.txt', delimiter=' ')
data = data.reshape((-1,10,12))
Lx = 82.4293
Lz = 81.004

plt.rc('text', usetex=True)
plt.rc('font', family='serif', weight='bold')
plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
plt.rc('font', size=28)
plt.rc('xtick', labelsize=28)
plt.rc('ytick', labelsize=28)

data_txz = np.zeros((data.shape[0], data.shape[1], data.shape[2]))
data_txz[:, :, :] = np.copy(data[:, :, :])
ligdata = data_txz
T = len(ligdata)

# get theta values at each lattice site as a function of time
theta_lat = []
for i in range(T):
    theta_lat.append( ligdata[i, :, :].flatten() )
theta_lat = np.array(theta_lat)
theta_lat = theta_lat.reshape((-1, 10, 12))
print np.mean(theta_lat)
# print theta_lat.shape

# plt.clf();plt.hist(theta_lat[:,:].flatten(), bins=100, normed=True);plt.show()
# print theta_lat[:,:].flatten().var()

theta_mean = np.mean(theta_lat, axis=0)
# print theta_mean.shape
# print 'x-mean', np.mean(theta_mean, axis=0)
# print 'z-mean', np.mean(theta_mean, axis=1)
plt.imshow(theta_mean, aspect=1.2, cmap="seismic_r", origin="lower", interpolation="none", vmin=-0.4, vmax=-0.1)
# plt.imshow(theta_mean, aspect=1.2, cmap="seismic_r", origin="lower", interpolation="none", vmin=-0.8, vmax=-0.1)
plt.xticks(np.arange(0, 12, 1))
plt.yticks(np.arange(0, 10, 1))
plt.xlim(-0.5,11.5)
plt.ylim(-0.5,9.5)
for i in np.arange(-0.5,12,1.0):
    plt.vlines(i, -0.5, 9.5, linestyle='solid', linewidth=2)
for i in np.arange(-0.5,10,1.0):
    plt.hlines(i, -0.5, 11.5, linestyle='solid', linewidth=2)
plt.colorbar()
plt.show()

theta_lat -= theta_lat.mean()

if args.savefigs:
    name_arr = range(0, theta_lat.shape[0], 10)
    name_arr = np.array(name_arr)
    for j in range(len(name_arr)):
    # for j in range(0, 100, 20):
        matr = theta_lat[name_arr[j]].transpose()
        plt.imshow(matr, aspect=1.2, cmap="seismic_r", origin="lower", interpolation="none", vmin=-0.4, vmax=-0.1)
#         plt.imshow(matr, aspect=1.7, cmap="PRGn_r", origin="lower", interpolation="none", vmin=-0.8, vmax=-0.1)
        plt.yticks(np.arange(0, 12, 1))
        plt.xticks(np.arange(0, 10, 1))
        plt.ylim(-0.5,11.5)
        plt.xlim(-0.5,9.5)
        for i in np.arange(-0.5,12,1.0):
            plt.hlines(i, -0.5, 9.5, linestyle='solid', linewidth=2)
        for i in np.arange(-0.5,10,1.0):
            plt.vlines(i, -0.5, 11.5, linestyle='solid', linewidth=2)
        plt.colorbar()
        plt.savefig('lat-' + args.bias + '-{:05d}.png'.format(j))
#         plt.savefig('lat-0.7250-{:05d}.png'.format(j))
        plt.clf()


all_corr_xz = []
samplez = 10
# DT = t_steps / samplez
X = theta_lat.shape[1]
Z = theta_lat.shape[2]
DT = theta_lat.shape[0] / samplez
# fig = plt.figure()
# ax  = fig.add_subplot(111, projection='3d')
var_list = []
for sample in xrange(samplez):
    To = DT*sample
    Tf = DT*(sample+1)
    sub_txz = theta_lat[To:Tf, :, :]
    
    cov_xz = np.zeros((X/2, Z/2))
    for t in xrange(DT):
       for x in xrange(X/2):
            for z in xrange(Z/2):
#                 print cov_xz.shape
#                 print sub_txz[t,x,z].shape
#                 print sub_txz[t,x : x + X/2,z : z + Z/2].shape
                cov_xz += sub_txz[t, x, z] * sub_txz[t, x : x + X/2, z : z + Z/2]
    
    cov_xz /= (DT * X/2 * Z/2 )
    var_list.append(cov_xz[0,0])
    corr_xz = cov_xz / cov_xz[0,0]
#     corr_xz = cov_xz
    
    x = range(X/2)
    z = range(Z/2)
    xv, zv = np.meshgrid(x, z)
    
    all_corr_xz.append(corr_xz)

all_corr_xz = np.array(all_corr_xz)
m_corr_xz = np.mean(all_corr_xz, axis=0)
d_corr_xz = np.std (all_corr_xz, axis=0) / np.sqrt(samplez)
# print np.mean(np.array(var_list))
# print m_corr_xz[1:,0].shape, range(1,X/2)
 
# correlation plots with first data point removed
# plt.errorbar(range(1, X/2), m_corr_xz[1:,0], d_corr_xz[1:,0], c='b', label="X", linewidth=2)
plt.errorbar(range(0, X/2), m_corr_xz[:,0], d_corr_xz[:,0], c='b', label="X", linewidth=2)
plt.hlines(0, -1, X/2, linestyles="dashed")
plt.xlim([-1.0,X/2])
plt.ylim(-0.1,0.2)
# plt.errorbar(range(1, Z/2), m_corr_xz[0,1:], d_corr_xz[0,1:], c='g', label="Z", linewidth=2)
plt.errorbar(range(0, Z/2), m_corr_xz[0,:], d_corr_xz[0,:], c='g', label="Z", linewidth=2)
plt.hlines(0, -1, Z/2, linestyles="dashed")
plt.xlabel('x or z', fontsize=30)
plt.ylabel('G(x, z)', fontsize=30)
plt.legend(loc='upper right', fontsize=30)
plt.show()

exit(1)

