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
plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
plt.rc('font', size=28)
plt.rc('font', weight='bold')
plt.rc('xtick', labelsize=28)
plt.rc('ytick', labelsize=28)

temp = 370.0
data = np.genfromtxt('f_i-{}.txt'.format(args.temp))
theta_axis = data[:,0]
f_i = data[:,1]
kB = 1.3806503 * 6.0221415 / 4184.0
target_beta = 1 / (kB * args.temp)
Lx = 81
N = 240

# fit gaussians to free energy minima to get partition functions
def spring(x, x0, k):
    return 0.5 * k * (x - x0) ** 2

ord_spring = np.where((theta_axis >= -0.74) * (theta_axis <=-0.68))
x_ord = theta_axis[ord_spring]
y_ord = f_i[ord_spring]
y_ord_min = min(y_ord)
# print x_ord, y_ord
popt, pcov = curve_fit(spring, x_ord, y_ord - y_ord_min, p0=(-0.73, 2500.0))
mean_ord = popt[0]
k_ord = popt[1]
v_ord = (1 / k_ord) * N
fit_ord = spring(x_ord, mean_ord, k_ord)
# print popt[0], popt[1]

disord_spring = np.where((theta_axis >= -0.50) * (theta_axis <=-0.11))
x_disord = theta_axis[disord_spring]
y_disord = f_i[disord_spring]
y_disord_min = min(y_disord)
# print x_disord, y_disord
popt, pcov = curve_fit(spring, x_disord, y_disord - y_disord_min, p0=(-0.32, 200.0))
mean_disord = popt[0]
k_disord = popt[1]
v_disord = (1 / k_disord) * N
fit_disord = spring(x_disord, mean_disord, k_disord)
# print popt[0], popt[1]

# plot original free energy and fits
plt.plot(theta_axis, f_i, marker='o', color='#006d2c')
plt.plot(theta_axis, f_i, color='#006d2c', linewidth=4, alpha=0.5)
# plt.plot(x_ord, fit_ord + y_ord_min, 'bv')
# plt.plot(x_disord, fit_disord + y_disord_min, 'g^')
plt.xlabel(r'$\langle\theta_z\rangle$', fontsize=28)
plt.ylabel(r'$F(\langle\theta_z\rangle) \textrm{ (in kcal/mol-rad)}$', fontsize=28)
plt.show()

# calculate and print partition functions and full free energies
Q_ord = integrate.simps( np.exp(-(fit_ord + y_ord_min)), x_ord )
fullF_ord = - np.log(Q_ord) / target_beta
E_ord = integrate.simps(data[:,2][ord_spring], x_ord)
S_ord = integrate.simps(data[:,3][ord_spring], x_ord)

Q_disord = integrate.simps( np.exp(-(fit_disord + y_disord_min)), x_disord )
fullF_disord = - np.log(Q_disord) / target_beta
E_disord = integrate.simps(data[:,2][disord_spring], x_disord)
S_disord = integrate.simps(data[:,3][disord_spring], x_disord)

trans_guess = (E_disord - E_ord) / (S_disord - S_ord)

print "order parameter values: ord = {} disord = {}".format(mean_ord, mean_disord)
print "partition functions: ord = {} disord = {}".format(Q_ord, Q_disord)
print "total free energies: ord = {} disord = {}".format(fullF_ord, fullF_disord)
print "total energies: ord = {} disord = {}".format(E_ord, E_disord)
print "total entropies: ord = {} disord = {}".format(S_ord, S_disord)
print "estimated transition by setting delta F = 0: {} K".format(trans_guess)

# optimising to find surface tension
data = np.genfromtxt('f_i-combined-375.6.txt')
theta_axis = data[:,0]
f_i = data[:,1]
trans_temp = 375.6
trans_beta = 1/(kB * trans_temp)

def gauss(x, x0, v, n):
#     return np.exp( -(x-x0)**2 / (2*v/n) ) / np.sqrt(2*np.pi*v/n)
    return np.exp( -(x-x0)**2 / (2*v/n) )

# fit to gaussians to get mean and variance
ord_spring = np.where((theta_axis >= -0.75) * (theta_axis <=-0.68))
x_ord = theta_axis[ord_spring]
y_ord = f_i[ord_spring]
y_ord_min = min(y_ord)
popt, pcov = curve_fit(spring, x_ord, y_ord - y_ord_min, p0=(-0.73, 2500.0))
mean_ord = popt[0]
v_ord = (1 / popt[1]) * N

y_disord = f_i[disord_spring]
y_disord_min = min(y_disord)
popt, pcov = curve_fit(spring, x_disord, y_disord - y_disord_min, p0=(-0.32, 200.0))
mean_disord = popt[0]
v_disord = (1 / popt[1]) * N

print y_ord_min, y_disord_min

prob_i = np.exp(-f_i)
p1 = gauss(theta_axis, mean_ord, v_ord, N) * np.exp(-y_ord_min) 
p2 = gauss(theta_axis, mean_disord, v_disord, N) * np.exp(-y_disord_min) 
# p1 = gauss(theta_axis, mean_ord, v_ord, N)
# p2 = gauss(theta_axis, mean_disord, v_disord, N)

# diagnostic plot #1
plt.plot(theta_axis, p1, 'bv', label='ordered prob dist')
plt.plot(theta_axis, p2, 'g^', label='disordered prob dist')
plt.plot(theta_axis, prob_i, 'ro', label='total prob dist')
plt.legend(loc='best')
plt.show()

# diagnostic plot #2
plt.plot(theta_axis,-np.log(p1), 'bv', label=r'\boldsymbol{-} \textbf{log(} \boldsymbol{P_{ord}} \textbf{)}')
plt.plot(theta_axis, -np.log(p2), 'g^', label=r'\boldsymbol{-} \textbf{log(} \boldsymbol{P_{disord}} \textbf{)}')
plt.plot(theta_axis, -np.log(prob_i), 'ro')
# plt.plot(theta_axis, -np.log(prob_i), 'r', linewidth=4, alpha=0.4)
plt.xlabel(r'\boldsymbol{$\langle\theta_z\rangle}', fontsize=28)
plt.ylabel(r'$\textbf{log probability}$', fontsize=28)
# plt.title(r'\textbf{Free energy at coexistence}')
plt.ylabel(r'$\boldsymbol{\beta F}$', fontsize=28)
plt.legend(loc='best', markerscale=2)
plt.show()

# diagnostic plot #3
g = prob_i - (p1 + p2)
plt.plot(theta_axis, p1, 'bv', label='ordered prob dist')
plt.plot(theta_axis, p2, 'g^', label='disordered prob dist')
plt.plot(theta_axis, prob_i, 'ro', label='total prob dist')
plt.plot(theta_axis, g, 'ks', label='difference in prob dist')
plt.legend(loc='best')
plt.show()

# # diagnostic plot #4
# plt.plot(theta_axis,-np.log(prob_i - (p1 + p2)), 'bv', label='log difference in prob dist')
# plt.plot(theta_axis, -np.log(prob_i), 'ro', label='log total prob dist')
# plt.legend(loc='best')
# plt.show()

# define fitting functions here
def surf_term(f, sigma):
    return np.exp( -trans_beta * sigma * Lx * min(2*np.pi*f, 1) )

def integrand(f, args):
    x, sigma = args
    num = -N * ( -mean_ord * f + mean_disord * (f - 1) + x )**2 / ( 2 * (f * (v_ord - v_disord) + v_disord) )
    denom = np.sqrt(2 * np.pi) * np.sqrt(v_disord / N) * np.sqrt( v_ord / (f*N - N * f**2) )
    denom = denom * np.sqrt( f * N * ( 1/v_ord + f/( v_disord * (1-f) ) ) )
    return (np.exp(num) * surf_term(f, sigma)) / denom 

def curve(x, sigma):
    res = integrate.quad(integrand, 0, 1, [x,sigma])
    ord_term = gauss(x, mean_ord, v_ord, N) * np.exp(-y_ord_min)
    disord_term = gauss(x, mean_disord, v_disord, N) * np.exp(-y_disord_min)
#     ord_term = gauss(x, mean_ord, v_ord, N)
#     disord_term = gauss(x, mean_disord, v_disord, N)
    return -np.log( res[0] + ord_term + disord_term )

vcurve = np.vectorize(curve)

# plot what free energy fit looks like for different values of sigma
size = np.logspace(np.log10(0.5), -1, 10).shape[0]
color = iter(plt.cm.rainbow(np.linspace(0,1,size)))
for sigma_i in np.logspace(np.log10(0.5), -1, 10):
    c = next(color)
    test_data = vcurve(theta_axis, sigma_i)
    plt.plot(theta_axis, test_data, 'o', color=c, label=r'$\textbf{{fit for }} \boldsymbol{{\sigma}} \textbf{{ = {:1.3f} kcal/mol-nm}}$'.format(sigma_i / (10 * trans_beta)))
    plt.plot(theta_axis, test_data, color=c, linewidth=4, alpha=0.5)
plt.plot(theta_axis, f_i, 'b', linewidth=3, alpha=0.8, label=r'$\textbf{original } \boldsymbol{\beta F}$')
plt.xlabel(r'$\boldsymbol{\langle\theta_z\rangle}', fontsize=28)
plt.ylabel(r'$\boldsymbol{\beta F(\langle\theta_z\rangle)}$', fontsize=28)
plt.legend(loc='best', fontsize=20)
plt.show()

# fitting procedure for best sigma value
popt, pcov = curve_fit(vcurve, theta_axis, f_i, p0=0.17)
print "results of surface tension optimisation (in kBT/Angstrom):"
print popt[0]
print "results of surface tension optimisation (in kcal/mol-nm):"
print popt[0] / (10 * trans_beta)

sigma_opt = popt[0]
test_data = vcurve(theta_axis, sigma_opt)
plt.plot(theta_axis, test_data, 'ro', label=r'$\textbf{{ fit for }} \boldsymbol{{\sigma}} \textbf{{ = {:1.3f} kcal/mol-nm}}'.format(sigma_opt / (10 * trans_beta)))
plt.plot(theta_axis, test_data, 'r', linewidth=4, alpha=0.5)
plt.plot(theta_axis, f_i, 'b', linewidth=3, alpha=0.5, label=r'$\textbf{original } \boldsymbol{\beta F}$')
plt.xlabel(r'$\boldsymbol{\langle\theta_z\rangle}', fontsize=28)
plt.ylabel(r'$\boldsymbol{\beta F(\langle\theta_z\rangle)}$', fontsize=28)
plt.legend(loc='best', fontsize=20)
plt.show()



