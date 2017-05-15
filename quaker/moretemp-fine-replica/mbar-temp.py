#!/usr/bin/python

import numpy as np
from math import *
import pymbar # for MBAR analysis
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("-temp", type=float, help="target temp at which to get free energy curve")
parser.add_argument("-show_err", action='store_true', help="to show free energy curve with error bars or not")
args = parser.parse_args()

free_energies_filename = 'f_k.out'

T = 2001		# number of snapshots
nbins_per_angle = 250	# number of bins per angle dimension
target_temperature = args.temp

def read_file(filename):
    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()
    return lines

# Create numpy array of temperatures

temperature_k = [349.5, 350.0, 350.5, 351.0, 351.5, 352.0, 352.5, 353.0]
temperature_k = np.array(temperature_k)
# print temperature_k
K = len(temperature_k)

# Compute inverse temperatures

kB = 1.3806503 * 6.0221415 / 4184.0 # Boltzmann constant in kcal/mol/K
beta_k = 1/ (kB * temperature_k)

# Read potential energies

U_kn = np.zeros((K, T))		# U_kn[k,t] is the potential energy (in kcal/mol) for snapshot t of temperature index k
for k in range(K):
    lines = read_file("pot-new.{:3.1f}.txt".format(temperature_k[k]))
    num = len(lines) - 2001
    U_kn[k, :] = lines[num:]

U_kn = np.array(U_kn)
# U_kn = np.transpose(U_kn)
# print U_kn.shape

# Read theta_z values

theta_kn = np.zeros((K, T))	# contains order parameter values in radians for snapshot t of temperature index k
for k in range(K):
    theta_k = np.genfromtxt("theta{:3.1f}.txt".format(temperature_k[k]))
    theta_k = theta_k.reshape((-1, 240))
    theta_lines = np.mean(theta_k, axis=1)
    num = len(theta_lines) - 2001
    theta_kn[k, :] = theta_lines[num:]

# theta_kn = np.array(theta_kn)
# print theta_kn.shape

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
for i in range(nbins_per_angle):
      val = angle_min + dx * (i + 0.5)
      # Determine which configurations lie in this bin.
      in_bin = (val-dx/2 <= theta_kn[indices]) & (theta_kn[indices] < val+dx/2)

      # Count number of configurations in this bin.
      bin_count = in_bin.sum()

      # Generate list of indices in bin.
      indices_in_bin = (indices[0][in_bin], indices[1][in_bin])

      if (bin_count > 0):
         bin_centers.append( (val) )
         bin_counts.append( bin_count )

         # assign these conformations to the bin index
         bin_kn[indices_in_bin] = nbins

         # increment number of bins
         nbins += 1

# print bin_centers

# Initialize MBAR

mbar = pymbar.MBAR(u_kln, N_k, verbose=False)
Deltaf_ij, dDeltaf_ij, Theta_ij = mbar.getFreeEnergyDifferences()
# print Deltaf_ij

[Delta_f_ij, dDelta_f_ij, Delta_u_ij, dDelta_u_ij, Delta_s_ij, dDelta_s_ij] = mbar.computeEntropyAndEnthalpy()
print Delta_s_ij

# plt.plot(temperature_k, Delta_u_ij[0], 'o')
# plt.show()

# Compute PMF at the desired temperature.

target_beta = 1.0 / (kB * target_temperature)
u_kn = target_beta * U_kn
(f_i, df_i) = mbar.computePMF(u_kn, bin_kn, nbins, uncertainties='from-lowest')
# print f_i

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xlabel(r'$\langle\theta_z\rangle$', fontsize=32)
plt.ylabel(r'-\beta\Delta F(\langle\theta_z\rangle)$', fontsize=32)
plt.xticks(fontsize=32, fontweight='bold')
plt.yticks(fontsize=32, fontweight='bold')

if args.show_err:
    plt.fill_between(bin_centers, f_i - 2*df_i, f_i + 2*df_i, alpha=.4)
else:
    plt.plot(bin_centers, f_i, color="#2020CC", linewidth=4)
plt.show()
