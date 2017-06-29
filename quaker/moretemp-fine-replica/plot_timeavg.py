import numpy as np
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="")
parser.add_argument("-temp", type=str, help="temp value to analyse")
parser.add_argument("-step", type=int, default=500, help="step value for running mean")
args = parser.parse_args()

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / N

data = np.genfromtxt('theta' + args.temp + '.txt', delimiter=' ')
data = -1.0 * data
data = data.reshape((-1,20,12))

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

data_txz = np.zeros(data.shape)
data_txz[:, ::2, :] = data[:, 0:10, :]
data_txz[:, 1::2, :] = data[:, 10:20, :]
data = data_txz

mean_tz = np.mean(data, axis=1)
mean_t = np.mean(mean_tz, axis=1)
plt.figure()
plt.plot(running_mean(mean_t, args.step), 'ro')
plt.show()
