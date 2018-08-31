import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=28)
plt.rc('font', weight='bold')
plt.rc('xtick', labelsize=28)
plt.rc('ytick', labelsize=28)

T = 25000

mean = 0.0
data = np.random.random_sample(T) - 0.5
data += mean

time_arr = np.linspace(0, T, T)
ACF = np.correlate(data, data, mode='full')
tcorr_0 = ACF[ACF.size/2:]
tcorr_0 /= tcorr_0[0]
# print time_arr.shape, tcorr.shape

mean = 4.5
data = np.random.random_sample(T) - 0.5
data += mean

time_arr = np.linspace(0, T, T)
ACF = np.correlate(data, data, mode='full')
tcorr_5 = ACF[ACF.size/2:]
tcorr_5 /= tcorr_5[0]
# print time_arr.shape, tcorr.shape

plt.plot(time_arr, tcorr_0, 'ro', label='mean = 0.0')
plt.plot(time_arr, tcorr_0, 'r', linewidth=4, alpha=0.5)
plt.plot(time_arr, tcorr_5, 'bv', label='mean = 5.0')
plt.plot(time_arr, tcorr_5, 'b', linewidth=4, alpha=0.5)
plt.legend(loc='best')
plt.show()


