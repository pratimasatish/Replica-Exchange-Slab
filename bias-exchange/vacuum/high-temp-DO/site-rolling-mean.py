import numpy as np
import argparse
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math

parser = argparse.ArgumentParser(description="")
parser.add_argument("-bias", type=str, help="bias value to analyse")
parser.add_argument("-step", type=int, default=500, help="step value for running mean")
args = parser.parse_args()

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / N

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

time_arr = np.arange(T)
# print running_mean(theta_lat[:,2,4], args.step).shape
# exit(1)

plt.plot(running_mean(theta_lat[:,2,4], args.step), 'r')
plt.plot(running_mean(theta_lat[:,13,8], args.step), 'b')
plt.show()
plt.plot(running_mean(theta_lat[:,4,10], args.step), 'r')
plt.plot(running_mean(theta_lat[:,5,1], args.step), 'b')
plt.show()
plt.plot(running_mean(theta_lat[:,6,3], args.step), 'r')
plt.plot(running_mean(theta_lat[:,17,9], args.step), 'b')
plt.show()

