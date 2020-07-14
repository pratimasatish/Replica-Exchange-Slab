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

data = np.genfromtxt('theta' + args.bias + '.txt', delimiter=' ')
data = data.reshape((-1,20,12))
Lx = 82.4293
Lz = 81.004

plt.rc('text', usetex=True)
plt.rc('font', family='serif', weight='bold')
plt.rc('xtick', labelsize=24)
plt.rc('ytick', labelsize=24)

data_txz = np.zeros((data.shape[0], data.shape[1], data.shape[2]))
data_txz[:, ::2, :] = data[:, 0:10, :]
data_txz[:, 1::2, :] = data[:, 10:20, :]
ligdata = data_txz
T = len(ligdata)

# get theta values at each lattice site as a function of time
theta_lat = []
for i in range(T):
    theta_lat.append( ligdata[i, :, :].flatten() )
theta_lat = np.array(theta_lat)
theta_lat = theta_lat.reshape((-1, 20, 12))
print theta_lat.shape

# plt.clf();plt.hist(theta_lat[:,:].flatten(), bins=100, normed=True);plt.show()

name_arr = [ [2,4], [3,4], [4,4], [5,4], [6,4] ]
name_arr = np.array(name_arr)
time_arr = np.arange(T)

plt.plot(time_arr, theta_lat[:,2,4], 'r')
plt.plot(time_arr, theta_lat[:,3,4], 'b')
plt.show()
plt.plot(time_arr, theta_lat[:,4,4], 'r')
plt.plot(time_arr, theta_lat[:,5,4], 'b')
plt.show()
plt.plot(time_arr, theta_lat[:,6,4], 'r')
plt.plot(time_arr, theta_lat[:,7,4], 'b')
plt.show()

fig = plt.figure(tight_layout=True, figsize=(10,9))
fig.subplots_adjust(hspace=0.0, wspace=0.0)
subplot = fig.add_subplot(2, 1, 1)
subplot.plot(time_arr, theta_lat[:,2,4], color='r')
subplot = fig.add_subplot(2, 1, 2)
subplot.plot(time_arr, theta_lat[:,3,4], color='b')

plt.show()

subplot = fig.add_subplot(2, 1, 1)
subplot.plot(time_arr, theta_lat[:,4,4], color='k')
subplot = fig.add_subplot(2, 1, 2)
subplot.plot(time_arr, theta_lat[:,5,4], color='g')

plt.show()

subplot = fig.add_subplot(2, 1, 1)
subplot.plot(time_arr, theta_lat[:,6,4], color='c')
subplot = fig.add_subplot(2, 1, 2)
subplot.plot(time_arr, theta_lat[:,7,4], color='m')

plt.show()

# name_arr = np.zeros((theta_lat.shape[1]*theta_lat.shape[2], 2))
# count = 0
# for i in range(theta_lat.shape[1]):
#     for j in range(theta_lat.shape[2]):
#         name_arr[count][0] = i
#         name_arr[count][1] = j
#         count = count + 1
# 
# color=iter(plt.cm.rainbow(np.linspace(0,1,name_arr.shape[0])))
# for i in name_arr:
#     c = next(color)
# #     plt.clf()
#     plt.plot(time_arr, theta_lat[:,int(i[0]), int(i[1])], color=c)
#     plt.title('site={},{}'.format(int(i[0]), int(i[1])))
#     plt.show()
# 
# exit(1)

theta_lat = []
for i in range(T):
    theta_lat.append( ligdata[i, :, :].flatten() )
theta_lat = np.array(theta_lat)


