#!/usr/bin/python

import numpy as np
from math import *
import pymbar # for MBAR analysis
from pymbar import timeseries # for timeseries analysis
import commands
import os
import os.path

free_energies_filename = 'f_k.out'

T = 2001		# number of snapshots
nbins_per_torsion = 10	# number of bins per torsion dimension
target_temperature = 351.0

def read_file(filename):
    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()
    return lines

# Create numpy array of temperatures

temperature_k = [349.0, 349.5, 350.0, 350.5, 351.0, 351.5]
temperature_k = np.array(temperature_k)
print temperature_k
K = len(temperature_k)

# Compute inverse temperatures

kB = 1.3806503 * 6.0221415 / 4184.0 # Boltzmann constant in kcal/mol/K
beta_k = 1/ (kB * temperature_k)

# Read potential eneriges

U_kn = np.zeros((K, T))		# U_kn[k,t] is the potential energy (in kcal/mol) for snapshot t of temperature index k
for k in range(K):
    lines = read_file("pot-new.{:3.1f}.txt".format(temperature_k[k]))
    U_kn[k, :] = lines[:]

U_kn = np.array(U_kn)
# U_kn = np.transpose(U_kn)
print U_kn.shape

# Read theta_z values

theta_kn = np.zeros((K, T))	# contains order parameter values in radians for snapshot t of temperature index k
for k in range(K):
    theta_k = np.genfromtxt("theta{:3.1f}.txt".format(temperature_k[k]))
    theta_k = theta_k.reshape((-1, 240))
    theta_lines = np.mean(theta_k, axis=1)
    theta_kn[k, :] = theta_lines[:]

# theta_kn = np.array(theta_kn)
print theta_kn.shape

# N_k = np.zeros([K], np.int32)
# N_k[:] = T
# N_max = T
# 
# # Create a list of indices of all configurations in kn-indexing.
# 
# mask_kn = np.zeros([K,N_max], dtype=np.bool)
# for k in range(0,K):
#    mask_kn[k,0:N_k[k]] = True
# # Create a list from this mask.
# indices = np.where(mask_kn)
# 
# # Compute reduced potential energy of all snapshots at all temperatures
# 
# u_kln = np.zeros([K,K,N_max], np.float32)		# u_kln[k,l,n] is reduced potential energy of trajectory segment n of temperature k evaluated at temperature l
# for k in range(K):
#    for l in range(K):
#       u_kln[k,l,0:N_k[k]] = beta_k[l] * U_kn[k,0:N_k[k]]
# 
# # Bin torsions into histogram bins for PMF calculation
# 
# angle_min = -np.pi
# angle_max = +np.pi
# dx = (angle_max - angle_min) / float(nbins_per_angle)
# # Assign angle bins
# bin_kn = np.zeros([K,N_max], np.int16)		# bin_kn[k,n] is the index of which histogram bin sample n from temperature index k belongs to
# nbins = 0
# bin_counts = list()
# bin_centers = list()		# bin_centers[i] is a theta_z value that gives the center of bin i
# # indices = np.arange(T)
# for i in range(nbins_per_torsion):
#       val = angle_min + dx * (i + 0.5)
#       # Determine which configurations lie in this bin.
#       in_bin = (val-dx/2 <= theta_kn[indices]) & (theta_kn[indices] < val+dx/2)
# 
#       # Count number of configurations in this bin.
#       bin_count = in_bin.sum()
# 
#       # Generate list of indices in bin.
#       indices_in_bin = (indices[in_bin])
# 
#       if (bin_count > 0):
#          bin_centers.append( (val) )
#          bin_counts.append( bin_count )
# 
#          # assign these conformations to the bin index
#          bin_kn[indices_in_bin] = nbins
# 
#          # increment number of bins
#          nbins += 1
# 
# # Initialize MBAR
# 
# mbar = pymbar.MBAR(u_kln, N_k, verbose=True)
# 
# # Compute PMF at the desired temperature.
# 
# target_beta = 1.0 / (kB * target_temperature)
# u_kn = target_beta * U_kn
# (f_i, df_i) = mbar.computePMF(u_kn, bin_kn, nbins, uncertainties='from-lowest')

