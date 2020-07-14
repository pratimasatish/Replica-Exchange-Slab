import numpy as np
import argparse
import scipy.integrate as integrate
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy import interpolate
import scipy as sp

parser = argparse.ArgumentParser(description="")
parser.add_argument("-temp", type=int, default=365, help="temperature")
args = parser.parse_args()

lwidth = 4.0
plt.rc('text', usetex=True, fontsize=26)
rcParams['text.latex.preamble'] = [r'\usepackage{helvet} \usepackage{sfmath}', r'\usepackage{upgreek}' ]
rcParams['axes.linewidth'] = 1.0*lwidth
rcParams['xtick.major.width'] = 1.0*lwidth
rcParams['xtick.major.size']  = 2.0*lwidth
rcParams['ytick.major.width'] = 1.0*lwidth
rcParams['ytick.major.size']  = 2.0*lwidth
plt.rc('lines', linewidth=4)
plt.rc('legend', frameon=False)  

temp = 365.0
data = np.genfromtxt('f_i-{}.txt'.format(args.temp))
theta_axis = data[:,0]
f_i = data[:,1]
prob_i = np.exp(-f_i)
kB = 1.3806503 * 6.0221415 / 4184.0
target_beta = 1 / (kB * args.temp)
Lx = 8.1
N = 240

# fit gaussians to free energy minima to get partition functions
def gaussian(x, x0, v, a):
    return a * np.exp(- (0.5 * (x - x0)**2) / v )

def parabola(x, x0, k, a, b):
    return a * 0.5 * k *(x - x0)**2 + b

ord_spring = np.where((theta_axis >= -0.76) * (theta_axis <=-0.68))
x_ord = theta_axis[ord_spring]
disord_spring = np.where((theta_axis >= -0.50) * (theta_axis <=-0.11))
x_disord = theta_axis[disord_spring]

# # gaussian fit to probability
# y_ord = prob_i[ord_spring]
# # print x_ord, y_ord
# popt, pcov = curve_fit( gaussian, x_ord, y_ord, p0=( np.mean(y_ord), np.var(y_ord), max(y_ord) ) )
# mean_ord = popt[0]
# v_ord = popt[1]
# norm_ord = popt[2]
# fit_ord = gaussian(x_ord, mean_ord, v_ord, norm_ord)
# print popt[0], popt[1], popt[2]
# 
# y_disord = prob_i[disord_spring]
# # print x_disord, y_disord
# popt, pcov = curve_fit(gaussian, x_disord, y_disord, p0=( np.mean(y_disord), np.var(y_disord), max(y_disord) ) )
# mean_disord = popt[0]
# v_disord = popt[1]
# norm_disord = popt[2]
# fit_disord = gaussian(x_disord, mean_disord, v_disord, norm_disord)
# print popt[0], popt[1], popt[2]

# spring fit to free energy
y_ord = f_i[ord_spring]
# print x_ord, y_ord
popt, pcov = curve_fit( parabola, x_ord, y_ord, p0=( np.mean(x_ord), 1/np.var(y_ord), 10, 0 ) )
mean_ord = popt[0]
k_ord = popt[1]
norm_ord = popt[2]
shift_ord = popt[3]
fit_ord = parabola(x_ord, mean_ord, k_ord, norm_ord, shift_ord)
print popt[0], popt[1], popt[2]

y_disord = f_i[disord_spring]
# print x_disord, y_disord
popt, pcov = curve_fit( parabola, x_disord, y_disord, p0=( np.mean(x_disord), 1/np.var(y_disord), 10, 10 ) )
mean_disord = popt[0]
k_disord = popt[1]
norm_disord = popt[2]
shift_disord = popt[3]
fit_disord = parabola(x_disord, mean_disord, k_disord, norm_disord, shift_disord)
print popt[0], popt[1], popt[2]

# plot original free energy and fits

# # gaussian fits
# plt.plot(theta_axis, prob_i, marker='o', color='#8856a7')
# plt.plot(theta_axis, prob_i, color='#8856a7', linewidth=4, alpha=0.5)

# spring fits
plt.plot(theta_axis, f_i, marker='o', color='#8856a7')
plt.plot(theta_axis, f_i, color='#8856a7', linewidth=4, alpha=0.5)
plt.plot(x_ord, fit_ord, 'bv')
plt.plot(x_disord, fit_disord, 'g^')
plt.xlabel(r'$\langle\theta_z\rangle$')
plt.ylabel(r'$F(\langle\theta_z\rangle) \textrm{ (in kcal/mol-rad)}$')
plt.show()

# calculate and print partition functions and full free energies

# for gaussian fits
Q_ord = integrate.simps( fit_ord, x_ord )
fullF_ord = - np.log(Q_ord) / target_beta
Q_disord = integrate.simps( fit_disord, x_disord )
fullF_disord = - np.log(Q_disord) / target_beta

# for spring fits
Q_ord = integrate.simps( np.exp(-fit_ord * target_beta), x_ord )
fullF_ord = - np.log(Q_ord) / target_beta
Q_disord = integrate.simps( np.exp(-fit_disord * target_beta), x_disord )
fullF_disord = - np.log(Q_disord) / target_beta

E_ord = integrate.simps(data[:,2][ord_spring], x_ord)
S_ord = integrate.simps(data[:,3][ord_spring], x_ord)

E_disord = integrate.simps(data[:,2][disord_spring], x_disord)
S_disord = integrate.simps(data[:,3][disord_spring], x_disord)

trans_guess = (E_disord - E_ord) / (S_disord - S_ord)
delS_acc = (E_disord - E_ord) / args.temp + kB * np.log(Q_disord / Q_ord)

print "order parameter values: ord = {} disord = {}".format(mean_ord, mean_disord)
print "partition functions: ord = {} disord = {}".format(Q_ord, Q_disord)
print "total free energies: ord = {} disord = {}".format(fullF_ord, fullF_disord)
print "total energies: ord = {} disord = {}".format(E_ord, E_disord)
print "total entropies: ord = {} disord = {}".format(S_ord, S_disord)
print "estimated transition by setting delta F = 0: {} K".format(trans_guess)
print "differences: delta F = {} delta E = {} delta S = {}".format(fullF_disord-fullF_ord, E_disord-E_ord, S_disord-S_ord)
print "entropy difference with more accurate method: {}".format(delS_acc)

# optimising to find surface tension
data = np.genfromtxt('f_i-368-from-365.txt')
theta_axis = data[:,0]
f_i = data[:,1]
prob_i = np.exp(-f_i)
trans_temp = 375.6
trans_beta = 1/(kB * trans_temp)

# fit to gaussians to get mean and variance
ord_spring = np.where((theta_axis >= -0.75) * (theta_axis <=-0.68))
x_ord = theta_axis[ord_spring]
y_ord = prob_i[ord_spring]
popt, pcov = curve_fit(gaussian, x_ord, y_ord, p0=( np.mean(y_ord), np.var(y_ord), max(y_ord) ))
mean_ord = popt[0]
v_ord = popt[1]
v_fit_ord = v_ord * N
norm_ord = popt[2]

disord_spring = np.where((theta_axis >= -0.50) * (theta_axis <=-0.11))
x_disord = theta_axis[disord_spring]
y_disord = prob_i[disord_spring]
popt, pcov = curve_fit(gaussian, x_disord, y_disord, p0=( np.mean(y_ord), np.var(y_ord), max(y_ord) ))
mean_disord = popt[0]
v_disord = popt[1]
v_fit_disord = v_disord * N
norm_disord = popt[2]

print "transition temp minima: ord = {} disord = {}".format(mean_ord, mean_disord)
print "transition temp sigma: ord = {} disord = {}".format(v_ord, v_disord)

p1 = gaussian(theta_axis, mean_ord, v_ord, norm_ord)
p2 = gaussian(theta_axis, mean_disord, v_disord, norm_disord)

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
    return np.exp( -trans_beta * sigma * Lx )
#     if (f <= 0.5):
#         return np.exp( -trans_beta * sigma * Lx * min(2*np.pi*f, 1) )
#     else:
#         return np.exp( -trans_beta * sigma * Lx * min(2*np.pi*(1-f), 1) )

def integrand(f, args):
    x, sigma = args
    num = -N * ( -mean_ord * f + mean_disord * (f - 1) + x )**2 / ( 2 * (f * (v_fit_ord - v_fit_disord) + v_fit_disord) )
    denom = np.sqrt(2 * np.pi) * np.sqrt(v_fit_disord / N) * np.sqrt( v_fit_ord / (f*N - N * f**2) )
    denom = denom * np.sqrt( f * N * ( 1/v_fit_ord + f/( v_fit_disord * (1-f) ) ) )
    return (np.exp(num) * surf_term(f, sigma)) / denom 

def curve(x, sigma):
    res = integrate.quad(integrand, 0.15, 0.9, [x,sigma])
#     res = integrate.quad(integrand, 0, 1, [x,sigma])
    ord_term = gaussian(x, mean_ord, v_ord, norm_ord)
    disord_term = gaussian(x, mean_disord, v_disord, norm_disord)
#     return -np.log( res[0] )
    return -np.log( res[0] + ord_term + disord_term )

vcurve = np.vectorize(curve)

# plot what free energy fit looks like for different values of sigma
size = np.logspace(np.log10(3.5), -0.5, 10).shape[0]
color = iter(plt.cm.rainbow(np.linspace(0,1,size)))
for sigma_i in np.logspace(np.log10(3.5), -0.5, 10):
    c = next(color)
    test_data = vcurve(theta_axis, sigma_i)
    plt.plot(theta_axis, test_data, 'o', color=c, label=r'$\textbf{{fit for }} \boldsymbol{{\sigma}} \textbf{{ = {:1.3f} kcal/mol-nm}}$'.format( sigma_i / trans_beta ))
    plt.plot(theta_axis, test_data, color=c, linewidth=4, alpha=0.5)
plt.plot(theta_axis, f_i, 'b', linewidth=3, alpha=0.8, label=r'$\textbf{original } \boldsymbol{\beta F}$')
plt.plot(theta_axis, -np.log(gaussian(theta_axis, mean_ord, v_ord, norm_ord)), 'kv')
plt.plot(theta_axis, -np.log(gaussian(theta_axis, mean_ord, v_ord, norm_ord)), 'k', linewidth=3, alpha=0.7)
plt.plot(theta_axis, -np.log(gaussian(theta_axis, mean_disord, v_disord, norm_disord)), 'kv')
plt.plot(theta_axis, -np.log(gaussian(theta_axis, mean_disord, v_disord, norm_disord)), 'k', linewidth=3, alpha=0.7)
plt.xlabel(r'$\boldsymbol{\langle\theta_z\rangle}', fontsize=28)
plt.ylabel(r'$\boldsymbol{\beta F(\langle\theta_z\rangle)}$', fontsize=28)
plt.legend(loc='best', fontsize=20)
plt.show()

size = np.logspace(np.log10(3.5), -0.5, 10).shape[0]
color = iter(plt.cm.rainbow(np.linspace(0,1,size)))
for sigma_i in np.logspace(np.log10(3.5), -0.5, 10):
    c = next(color)
    test_data = vcurve(theta_axis, sigma_i)
    plt.plot(theta_axis, np.exp(-test_data), 'o', color=c, label=r'$\textbf{{fit for }} \boldsymbol{{\sigma}} \textbf{{ = {:1.3f} kcal/mol-nm}}$'.format( sigma_i / trans_beta ))
    plt.plot(theta_axis, np.exp(-test_data), color=c, linewidth=4, alpha=0.5)
plt.plot(theta_axis, prob_i, 'b', linewidth=3, alpha=0.8, label=r'$\textbf{original } \boldsymbol{P(\theta)}$')
# plt.plot(theta_axis, gaussian(theta_axis, mean_ord, v_ord, norm_ord), 'kv')
# plt.plot(theta_axis, gaussian(theta_axis, mean_ord, v_ord, norm_ord), 'k', linewidth=3, alpha=0.7)
# plt.plot(theta_axis, gaussian(theta_axis, mean_disord, v_disord, norm_disord), 'kv')
# plt.plot(theta_axis, gaussian(theta_axis, mean_disord, v_disord, norm_disord), 'k', linewidth=3, alpha=0.7)
plt.xlabel(r'$\boldsymbol{\langle\theta_z\rangle}')
plt.ylabel(r'$\boldsymbol{\beta F(\langle\theta_z\rangle)}$')
plt.legend(loc='best', fontsize=20)
plt.show()


# fitting procedure for best sigma value
popt, pcov = curve_fit(vcurve, theta_axis, f_i, p0=0.17)
print "results of surface tension optimisation (in kBT/nm):"
print popt[0]
print "results of surface tension optimisation (in kcal/mol-nm):"
print popt[0] / (trans_beta)

sigma_opt = popt[0]
test_data = vcurve(theta_axis, sigma_opt)
plt.figure(figsize=(14,9))
#plt.plot(theta_axis, test_data, 'r', linewidth=4, alpha=0.5)
plt.plot(theta_axis, f_i, 'k', linewidth=4, label=r'$\mathsf{free\,\,energy\,\,density\,\,at\,\,}T_t$')
plt.plot(theta_axis, test_data, 'r', marker='v', markersize=10, linewidth=4, alpha=0.5, label=r'$\mathsf{{fit\,\,for\,\,}}\sigma ={:1.3f}\mathsf{{\,\,kcal/mol}}\mbox{{-}}\mathsf{{nm}}$'.format( sigma_opt / trans_beta ))
plt.xlabel(r'$\Theta_z\,\mathsf{(rad)}$')
plt.ylabel(r'$\mathsf{free\,\,energy\,\,density\,\,(kcal/mol}\mbox{-}\mathsf{rad)}$')
plt.xlim(-0.8, 0.0)
plt.legend(loc='upper right')
plt.savefig('best-surf-tens-coex-solv.pdf')






