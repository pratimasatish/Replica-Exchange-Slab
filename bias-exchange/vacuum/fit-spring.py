import numpy as np
import argparse
import scipy.integrate as integrate
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy import interpolate
import scipy as sp

parser = argparse.ArgumentParser(description="")
parser.add_argument("-temp", type=int, default=370, help="temperature")
args = parser.parse_args()

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=28)
plt.rc('font', weight='bold')
plt.rc('xtick', labelsize=28)
plt.rc('ytick', labelsize=28)

data = np.genfromtxt('f_i-{}.txt'.format(args.temp))
theta_axis = data[:,0]
f_i = data[:,1]
kB = 1.3806503 * 6.0221415 / 4184.0
target_beta = 1 / (kB * args.temp)

# fit gaussians to free energy minima to get partition functions
def spring(x, x0, k):
    return 0.5 * k * (x - x0) ** 2

ord_spring = np.where((theta_axis >= -0.74) * (theta_axis <=-0.67))
x_ord = theta_axis[ord_spring]
y_ord = f_i[ord_spring]
y_ord_min = min(y_ord)
# print x_ord, y_ord
popt, pcov = curve_fit(spring, x_ord, y_ord - y_ord_min, p0=(-0.73, 2500.0))
fit_ord = spring(x_ord, popt[0], popt[1])
# print popt[0], popt[1]

disord_spring = np.where((theta_axis >= -0.46) * (theta_axis <=-0.16))
x_disord = theta_axis[disord_spring]
y_disord = f_i[disord_spring]
y_disord_min = min(y_disord)
# print x_disord, y_disord
popt, pcov = curve_fit(spring, x_disord, y_disord - y_disord_min, p0=(-0.32, 200.0))
fit_disord = spring(x_disord, popt[0], popt[1])
# print popt[0], popt[1]

# plot original free energy and fits
plt.plot(theta_axis, f_i, 'ro')
plt.plot(theta_axis, f_i, 'r', linewidth=4, alpha=0.5)
plt.plot(x_ord, fit_ord + y_ord_min, 'bv')
plt.plot(x_disord, fit_disord + y_disord_min, 'g^')
plt.show()

# calculate and print partition functions and full free energies
Q_ord = integrate.simps( np.exp(-target_beta * (fit_ord + y_ord_min)), x_ord )
fullF_ord = - np.log(Q_ord) / target_beta

Q_disord = integrate.simps( np.exp(-target_beta * (fit_disord + y_disord_min)), x_disord )
fullF_disord = - np.log(Q_disord) / target_beta

print "partition functions: ord = {} disord = {}".format(Q_ord, Q_disord)
print "total free energies: ord = {} disord = {}".format(fullF_ord, fullF_disord)



