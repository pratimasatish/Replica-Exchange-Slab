import numpy as np
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("-t", type=float, help="temperature for which to calculate correlations")
args = parser.parse_args()

data = np.genfromtxt("theta{:3.1f}.txt".format(args.t))
# print np.sum(np.square(data), axis=1)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

n_lig = 240
n_steps = len(data)
tdata = np.linspace(0, n_steps, 200)

thz_avg = np.mean(data)
print thz_avg, n_steps

plt.figure(1)
xdata_t = data.reshape((-1,n_lig), order='C')
print xdata_t.shape
xcorr = np.zeros(n_steps)

for col in xrange(n_lig):
    lig = xdata_t[:,col]
#     print lig.shape
 
    for i in xrange(len(lig)):
        for j in xrange(n_steps-i):
            xcorr[j] = xcorr[j] + (lig[i] * lig[i+j] - thz_avg * thz_avg)

for i in xrange(len(xcorr)):
    xcorr[i] = xcorr[i] / ((n_steps-i) * n_lig)

plt.plot(tdata[1:], xcorr[1:], color="#fdc086", linewidth=2)
plt.show()



