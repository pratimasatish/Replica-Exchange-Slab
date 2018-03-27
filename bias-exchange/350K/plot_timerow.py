import numpy as np
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="")
parser.add_argument("-bias", type=str, help="bias value to analyse")
args = parser.parse_args()

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / N

data = np.genfromtxt('theta' + args.bias + '.txt', delimiter=' ')
data = -1.0 * data
data = np.mean(data.reshape((-1, 240)), axis=1)

log = np.genfromtxt('log.stripped')
time_max = len(data) * 25000
time_arr = np.array(range(0, time_max, 25000))

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.figure()
# plt.plot(log[:,1], data[1:], 'ro')
plt.plot(log[:-1,1], data, 'b')
# 
# plt.plot(time_arr, data, 'ro')
# plt.plot(time_arr, data, 'b')
plt.show()

bins = np.linspace(-1.70, 1.70, 100)
hist, bins = np.histogram(data, bins = bins, density = True)
bin_centres = bins[1:] * 0.5 + bins[:-1] * 0.5
plt.plot(bin_centres, hist)
plt.show()

# for i in np.arange(0.4, 1.6, 0.3):
#     hist, bins = np.histogram(data, bins = bins, density = True)
#     bin_centres = bins[1:] * 0.5 + bins[:-1] * 0.5

