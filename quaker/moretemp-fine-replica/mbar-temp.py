#!/usr/bin/python

import numpy as np
from math import *
import pymbar # for MBAR analysis
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import argparse
import matplotlib.tri
import scipy.integrate as integrate
from scipy.optimize import curve_fit

parser = argparse.ArgumentParser(description="")
parser.add_argument("-temp", type=float, help="target temp at which to get free energy curve")
parser.add_argument("-dim", type=float, default=1, help="order parameter dimensionality [1d = th_z, 2d = th_z and th_x]")
parser.add_argument("-show_err", action='store_true', help="to show free energy curve with error bars or not")
parser.add_argument("-prob", action='store_true', help="to show free energy or probability distribution")
args = parser.parse_args()

free_energies_filename = 'f_k.out'

T = 2001		# number of snapshots
nbins_per_angle = 500	# number of bins per angle dimension
target_temperature = args.temp
namelist_1 = np.arange(-0.8500, -0.0240, 0.0250)
namelist_2 = np.arange(0.0, 0.1040, 0.0250)
namelist = np.concatenate((namelist_1, namelist_2))
k_list = np.ones(len(namelist)) * 15000.0
N_bias = len(namelist) * 2
T_bias = 5000
# bias_theta_kn = np.zeros((N_sims, T))
# bias_theta_x_kn = np.zeros((N_sims, T))

# target_temperature = 350.1

def read_file(filename):
    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()
    return lines

# Create numpy array of temperatures

temperature_k = [349.5, 350.0, 350.5, 351.0, 351.5, 352.0, 352.5, 353.0]
temperature_k = np.array(temperature_k)
temperature_k = np.concatenate( (np.ones(N_bias/2) * 350.18, np.ones(N_bias/2) * 355, temperature_k) )
# temperature_k = np.concatenate((temperature_k, np.ones(N_bias) * 355))
# print temperature_k
K = len(temperature_k)
U_kn = np.zeros((K, T))		# U_kn[k,t] is the potential energy (in kcal/mol) for snapshot t of temperature index k

# Compute inverse temperatures

kB = 1.3806503 * 6.0221415 / 4184.0 # Boltzmann constant in kcal/mol/K
beta_k = 1/ (kB * temperature_k)

# Read theta_z values

# first for biased simulations at 350.18K

theta_kn = np.zeros((K, T))	# contains order parameter values in radians for snapshot t of temperature index k
for k, biasval in enumerate(namelist):
    data = np.genfromtxt('theta-350.18.{:1.4f}.txt'.format(biasval))
    data = data.reshape((-1, 240))
    data_t = np.mean(data, axis=1)
    num = len(data_t) - 2001
    print k, biasval
    theta_kn[k, :] = data_t[num:]

# Next for biased simulations at 355K

for k, biasval in enumerate(namelist):
    data = np.genfromtxt('theta-355.{:1.4f}.txt'.format(biasval))
    data = data.reshape((-1, 240))
    data_t = np.mean(data, axis=1)
    num = len(data_t) - 2001
    print k+N_bias/2, biasval
    theta_kn[k + N_bias/2, :] = data_t[num:]

# Last for replica exchange data

for k in range(8):
    theta_k = np.genfromtxt("theta{:3.1f}.txt".format(temperature_k[k + N_bias]))
    theta_k = theta_k.reshape((-1, 240))
    theta_lines = np.mean(theta_k, axis=1)
    num = len(theta_lines) - 2001
    print k+N_bias
    theta_kn[k + N_bias, :] = theta_lines[num:]

# theta_kn = np.array(theta_kn)
# print theta_kn.shape

theta_x_kn = np.zeros((K, T))	# contains order parameter values in radians for snapshot t of temperature index k

# First for biased simulations at 350.18K

for k, biasval in enumerate(namelist):
    data_x = np.genfromtxt('theta_x-350.18.{:1.4f}.txt'.format(biasval))
    data_x = data_x.reshape((-1, 240))
    data_xt = np.mean(data_x, axis=1)
    num = len(data_xt) - 2001
    theta_x_kn[k, :] = data_xt[num:]

# Next for biased simulations at 355K

for k, biasval in enumerate(namelist):
    data_x = np.genfromtxt('theta_x-355.{:1.4f}.txt'.format(biasval))
    data_x = data_x.reshape((-1, 240))
    data_xt = np.mean(data_x, axis=1)
    num = len(data_xt) - 2001
    theta_x_kn[k + N_bias/2, :] = data_xt[num:]

# Last for replica exchange data

for k in range(8):
    theta_x_k = np.genfromtxt("theta_x{:3.1f}.txt".format(temperature_k[k + N_bias]))
    theta_x_k = theta_x_k.reshape((-1, 240))
    theta_x_lines = np.mean(theta_x_k, axis=1)
    num = len(theta_x_lines) - 2001
    theta_x_kn[k + N_bias, :] = theta_x_lines[num:]

# Read potential energies, first from umbrella sampling at 350.18K, then from replica exchange

# At 350.18K

for k, th in enumerate(namelist):
    lines = np.genfromtxt("pot-350.18.{:1.4f}".format(th))
    num = len(lines) - 2001
#     print k, th
    U_kn[k, :] = lines[num:] + 0.5 * k_list[k] * (theta_kn[k, :] - th) * (theta_kn[k, :] - th)

# At 355K

for k, th in enumerate(namelist):
    lines = np.genfromtxt("pot-355.{:1.4f}".format(th))
    num = len(lines) - 2001
#     print k, th
    U_kn[k + N_bias/2, :] = lines[num:] + 0.5 * k_list[k] * (theta_kn[k + N_bias/2, :] - th) * (theta_kn[k + N_bias/2, :] - th)

# From replica exchange

for k in range(8):
    lines = np.genfromtxt("pot-new.{:3.1f}.txt".format(temperature_k[k + N_bias]))
    num = len(lines) - 2001
    U_kn[k + N_bias, :] = lines[num:]

U_kn = np.array(U_kn)
# U_kn = np.transpose(U_kn)
# print U_kn.shape

# Initialise MBAR arrays

N_k = np.zeros([K], np.int32)
N_k[:] = T
N_max = T

# Create a list of indices of all configurations in kn-indexing.

mask_kn = np.zeros([K,N_max], dtype=np.bool)
for k in range(0,K):
   mask_kn[k,0:N_k[k]] = True
# Create a list from this mask.
indices = np.where(mask_kn)

# print indices

# Compute reduced potential energy of all snapshots at all temperatures

u_kln = np.zeros([K,K,N_max], np.float32)		# u_kln[k,l,n] is reduced potential energy of trajectory segment n of temperature k evaluated at temperature l
for k in range(K):
   for l in range(K):
      u_kln[k,l,0:N_k[k]] = beta_k[l] * U_kn[k,0:N_k[k]]

# print u_kln.shape

# Bin angles into histogram bins for PMF calculation

angle_min = -np.pi
angle_max = +np.pi
dx = (angle_max - angle_min) / float(nbins_per_angle)
# Assign angle bins
bin_kn = np.zeros([K,N_max], np.int16)		# bin_kn[k,n] is the index of which histogram bin sample n from temperature index k belongs to
nbins = 0
bin_counts = list()
bin_centers = list()		# bin_centers[i] is a theta_z value that gives the center of bin i
# indices = np.arange(T)

if (args.dim == 2):
    # Get bins for 2d order parameter space (theta_z and theta_x)
    for i in range(nbins_per_angle):
        for j in range(nbins_per_angle):
            val = angle_min + dx * (i + 0.5)
            val_x = angle_min + dx * (j + 0.5)
            # Determine which configurations lie in this bin.
            in_bin = (val-dx/2 <= theta_kn[indices]) & (theta_kn[indices] < val+dx/2) & (val_x-dx/2 <= theta_x_kn[indices]) & (theta_x_kn[indices] < val_x+dx/2)
      
            # Count number of configurations in this bin.
            bin_count = in_bin.sum()
      
            # Generate list of indices in bin.
            indices_in_bin = (indices[0][in_bin], indices[1][in_bin])
      
            if (bin_count > 0):
                bin_centers.append( (val, val_x) )
                bin_counts.append( bin_count )
      
                # assign these conformations to the bin index
                bin_kn[indices_in_bin] = nbins
      
                # increment number of bins
                nbins += 1

# print bin_centers

else:
    # Get bins for 1d order parameter space (theta_z and theta_x)
    for i in range(nbins_per_angle):
        val = angle_min + dx * (i + 0.5)
        # Determine which configurations lie in this bin.
        in_bin = (val-dx/2 <= theta_kn[indices]) & (theta_kn[indices] < val+dx/2) 
      
        # Count number of configurations in this bin.
        bin_count = in_bin.sum()
    
        # Generate list of indices in bin.
        indices_in_bin = (indices[0][in_bin], indices[1][in_bin])
    
        if (bin_count > 0):
            bin_centers.append( val )
            bin_counts.append( bin_count )
     
            # assign these conformations to the bin index
            bin_kn[indices_in_bin] = nbins
     
            # increment number of bins
            nbins += 1

# Initialize MBAR

mbar = pymbar.MBAR(u_kln, N_k, verbose=False)
Deltaf_ij, dDeltaf_ij, Theta_ij = mbar.getFreeEnergyDifferences()
# print Deltaf_ij

# [Delta_f_ij, dDelta_f_ij, Delta_u_ij, dDelta_u_ij, Delta_s_ij, dDelta_s_ij] = mbar.computeEntropyAndEnthalpy()
# print Delta_s_ij

# plt.plot(temperature_k, Delta_u_ij[0], 'o')
# plt.show()

# Compute PMF at the desired temperature.

target_beta = 1.0 / (kB * target_temperature)
u_kn = target_beta * U_kn
(f_i, df_i) = mbar.computePMF(u_kn, bin_kn, nbins, uncertainties='from-lowest')
print f_i, bin_centers
prob_i = np.exp(-f_i)
dprob_i = np.exp(-f_i) * np.std(f_i)
# print len(prob_i)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xlabel(r'$\langle\theta_z\rangle$', fontsize=32)
plt.ylabel(r'-\beta\Delta F(\langle\theta_z\rangle)$', fontsize=32)
plt.xticks(fontsize=32, fontweight='bold')
plt.yticks(fontsize=32, fontweight='bold')

def gauss(x, x0, s, a):
    return a * np.exp(-(x - x0)**2 / (2 * s * s))
# 
# x_ord = np.array(bin_centers[:4])
# y = np.array(prob_i[:4])
# popt, pcov = curve_fit(gauss, x_ord, y, p0=(-0.75, 0.25, 0.8))
# fit_ord = gauss(bin_centers[:4], popt[0], popt[1], popt[2])
# integrand = lambda x: popt[2] * np.exp(- (x - popt[0])**2 / 2 * popt[1]**2 )
# area_ord = integrate.quad(integrand, -0.80, -0.67)
# print popt
# 
# x_disord = np.array(bin_centers[4:])
# y = np.array(prob_i[4:])
# popt, pcov = curve_fit(gauss, x_disord, y, p0=(-0.25, 0.45, 0.8))
# fit_disord= gauss(bin_centers[4:], popt[0], popt[1], popt[2])
# integrand = lambda x: popt[2] * np.exp(- (x - popt[0])**2 / 2 * popt[1]**2 )
# area_disord = integrate.quad(integrand, -0.45, 0.0)
# print popt
# 
# print area_ord, area_disord
# 

if (args.dim == 2):
    fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    ax = fig.gca(projection='3d')
    bin_centers = np.array(bin_centers)
    # print bin_centers[:,0]
    # print len(bin_centers[:,0]), len(f_i)
    thz = bin_centers[:,0]
    thx = bin_centers[:,1]
    thZ, thX = np.meshgrid(thz, thx)
    cols = np.unique(thz).shape[0]
    print cols
    # print thZ.shape, thz.shape, f_i.shape
    print thZ.shape, thX.shape, f_i.shape

    # if args.prob:
    #     ax.plot_trisurf(thz, thx, f_i, antialiased=True, linewidth=0)

# if args.show_err:
#     if args.prob:
#         plt.fill_between(bin_centers, prob_i - 2*dprob_i, prob_i + 2*dprob_i, alpha=.4)
#     else:
#         plt.fill_between(bin_centers, f_i - 2*df_i, f_i + 2*df_i, alpha=.4)
# else:
#     if args.prob:
#         ax.plot_trisurf(thz, thx, f_i, antialiased=True, linewidth=0)
# #         plt.plot(bin_centers, prob_i, color="#2020CC", linewidth=4, alpha=.4)
# #         plt.plot(x_ord, fit_ord, color='r')
# #         plt.plot(x_disord, fit_disord, color='r')
# #     else:
# #         plt.plot(bin_centers, f_i, color="#2020CC", linewidth=4)
plt.show()

area_ord = integrate.simps(prob_i[:7], bin_centers[:7])
area_disord = integrate.simps(prob_i[7:], bin_centers[7:])

print area_ord, area_disord
