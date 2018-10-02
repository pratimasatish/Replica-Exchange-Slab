import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.optimize import curve_fit

kB = 1.3806503 * 6.0221415 / 4184.0
target_temp = 370
target_beta = 1/(kB*target_temp)
lx = 81
N = 240
a1 = -0.7218
v1 = 0.01359
a2 = -0.2844
v2 = 0.6769

data = np.genfromtxt('g_th.dat')
g_th = data[:,1]

data = np.genfromtxt('f_v_thz-351K.dat')
thz = data[:,0]
f_i = data[:,1]
prob_i = np.exp(-f_i)

def gauss(x, x0, v):
    return np.exp(-(x-x0)**2 / 2*(v/N)) / np.sqrt(2*np.pi*v/N)

def surf_term(f, gamma):
    return np.exp( -target_beta * gamma * lx * min(2*np.pi*f, 1) )

def integrand(f, args):
    x, gamma = args
    num = -N * ( -a1 * f + a2 * (f - 1) + x )**2 / ( 2 * (f * (v1 - v2) + v2) )
    denom = 2 * np.pi * np.sqrt(v2 / N) * np.sqrt( v1 / (f*N - N * f**2) )
    denom = denom * np.sqrt( ( (N * v1 - N * v2) * f**2 + f*N*v2 ) / ( 2 * np.pi * v1 * v2 * (1-f) ) )
    return (np.exp(num) * surf_term(f, gamma)) / denom 

#     exp_num = 34.6116*f*f*f + (158.261*x+10.3966)*f*f + f*(180.911*x*x - 55.3587*x - 30.3767) - 180.911*x*x - 102.902*x - 14.6327
#     exp_denom = (f-1)*(f-1.02049)
#     exp_term = np.exp(exp_num / exp_denom)
#     pf1 = -7.58853/(f+0.0280367*x-1)
#     pf2 = np.sqrt((1-f) / (f*f - 2.02049*f + 1.02049))
#     pf3 = f*f + f*(0.0280367*x - 2) - 0.0280367*x + 1
#     return pf1 * pf2 * pf3 * exp_term * np.exp( -target_beta * gamma * lx * min(2*np.pi*f, 1) )

p1 = gauss(thz, a1, v1)
p2 = gauss(thz, a2, v2)

def curve(x, gamma):
     res = integrate.quad(integrand, 0, 1, [x,gamma])
     return -np.log( res[0] + 0.02 * ( gauss(x, a1, v1) + gauss(x, a2, v2) ) )

vcurve = np.vectorize(curve)
popt, pcov = curve_fit(vcurve, thz, f_i, p0=1e-5)

plt.plot(thz, vcurve(thz, popt[0]), 'bs')
plt.plot(thz, f_i, 'ro')
plt.show()


y_guess = prob_i - 0.02 * (p1 + p2)






