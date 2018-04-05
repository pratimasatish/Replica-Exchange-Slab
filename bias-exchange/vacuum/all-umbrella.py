import pymbar
import argparse
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from matplotlib import cm
from scipy.optimize import curve_fit
from scipy import interpolate
import scipy as sp
from mpl_toolkits.mplot3d import Axes3D

# parser = argparse.ArgumentParser(description="")
# parser.add_argument("-dim", type=float, default=1, help="order parameter dimensionality [1d = th_z, 2d = th_z and th_x]")
# args = parser.parse_args()

# initialise list of temperatures 
kB = 1.3806503 * 6.0221415 / 4184.0

namelist_350 = np.array( [-0.750, -0.700, -0.675, -0.625, -0.600, -0.550, -0.500, -0.450, -0.425, -0.400, -0.375, -0.350, -0.275, -0.225, -0.175, -0.125] )
namelist_365 = np.array( [-0.750, -0.700, -0.675, -0.625, -0.600, -0.550, -0.500, -0.450, -0.425, -0.400, -0.375, -0.350, -0.275, -0.225, -0.175, -0.125] )
namelist_370 = np.array( [-0.750, -0.700, -0.650, -0.625, -0.600, -0.580, -0.560, -0.540, -0.520, -0.500, -0.425, -0.350, -0.275, -0.225, -0.175, -0.125] )
namelist_375 = np.array( [-0.750, -0.700, -0.650, -0.625, -0.600, -0.580, -0.560, -0.540, -0.520, -0.500, -0.425, -0.350, -0.275, -0.225, -0.175, -0.125] )

full_namelist = np.concatenate(( namelist_350, namelist_365, namelist_370, namelist_375 ))
T_bias = 5000
N_bias = len(namelist_370)

temp_list = np.concatenate(( np.ones(N_bias) * 350.0, np.ones(N_bias) * 365.0, np.ones(N_bias) * 370.0, np.ones(N_bias) * 375.0 ))
beta_list = 1/(kB * temp_list)
k_list = np.ones(len(full_namelist)) * 1500.0

N_sims = len(temp_list)
theta_ik = []
VO_ik = []
UO_ik = []
 
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=28)
plt.rc('font', weight='bold')
plt.rc('xtick', labelsize=28)
plt.rc('ytick', labelsize=28)

# build up theta matrix

for biasval in namelist_350:
    data = np.genfromtxt('350K/theta{:1.3f}.txt'.format(biasval))
    data = data.reshape((-1, 240))
    data_i = np.mean(data, axis=1)
    theta_ik.append( data_i[1:] )

for biasval in namelist_365:
    data = np.genfromtxt('365K/theta{:1.3f}.txt'.format(biasval))
    data = data.reshape((-1, 240))
    data_i = np.mean(data, axis=1)
    theta_ik.append( data_i[1:] )

for biasval in namelist_370:
    data = np.genfromtxt('370K/theta{:1.3f}.txt'.format(biasval))
    data = data.reshape((-1, 240))
    data_i = np.mean(data, axis=1)
    theta_ik.append( data_i )

for biasval in namelist_375:
    data = np.genfromtxt('375K/theta{:1.3f}.txt'.format(biasval))
    data = data.reshape((-1, 240))
    data_i = np.mean(data, axis=1)
    theta_ik.append( data_i )

# build up potential matrix

for k, th in enumerate(namelist_350):
    lines = np.genfromtxt("350K/pot-new.{:1.3f}".format(th))
    VO_ik.append( lines[5001:] )
    dtheta_i = np.array(theta_ik[k]) - th
    UO_ik.append( lines[5001:] - 0.5 * k_list[k] * np.square(dtheta_i) )

for k, th in enumerate(namelist_365):
    lines = np.genfromtxt("365K/pot-new.{:1.3f}".format(th))
    VO_ik.append( lines[5001:] )
    dtheta_i = np.array(theta_ik[k + N_bias]) - th
    UO_ik.append( lines[5001:] - 0.5 * k_list[k + N_bias] * np.square(dtheta_i) )

for k, th in enumerate(namelist_370):
    lines = np.genfromtxt("370K/pot-new.{:1.3f}".format(th))
    VO_ik.append( lines[5001:] )
    dtheta_i = np.array(theta_ik[k + 2*N_bias]) - th
    UO_ik.append( lines[5001:] - 0.5 * k_list[k + 2*N_bias] * np.square(dtheta_i) )

for k, th in enumerate(namelist_375):
    lines = np.genfromtxt("375K/pot-new.{:1.3f}".format(th))
    VO_ik.append( lines[5001:] )
    dtheta_i = np.array(theta_ik[k + 3*N_bias]) - th
    UO_ik.append( lines[5001:] - 0.5 * k_list[k + 3*N_bias] * np.square(dtheta_i) )

# calculate sizes for MBAR arrays
N_k = [ len(VO_i) for VO_i in VO_ik ]
N_k = np.array(N_k)
u_mbar = np.zeros((len(VO_ik), sum(N_k)))
K = u_mbar.shape[0]
N = u_mbar.shape[1]

# make numpy arrays from data
N_max = max(N_k)
th_ik = np.zeros([K, N_max])
uo_ik = np.zeros([K, N_max])
k = 0
for line1, line2 in zip(theta_ik, UO_ik):
    th_ik[k,0:len(line1)] = np.array(line1)
    uo_ik[k,0:len(line2)] = np.array(line2)
    k = k + 1

# go row by row to evaluate configuration energy at each temperature
for k in range(K):
    # populate off-diagonal blocks in MBAR array; go column by column, i.e. config by config
    for i in range(N_k[k]):
        dtheta = th_ik[k, i] - full_namelist
        print k, i
        u_mbar[ :, sum(N_k[:k]) + i ] = beta_list * ( uo_ik[k,i] + 0.5 * k_list * np.square(dtheta) )

my_mbar = pymbar.MBAR(u_mbar, N_k)

theta_n = [val for row in theta_ik for val in row]
theta_n = np.array(theta_n)

mask_kn = np.zeros([K,N_max], dtype=np.bool)
for k in range(0,K):
   mask_kn[k,0:N_k[k]] = True
# Create a list from this mask.
indices = np.where(mask_kn)
max_bins = 100
nbins = 0
min_val = theta_n.min()
max_val = theta_n.max()
delta = (max_val - min_val) / float(max_bins)
bin_kn = np.zeros([K,N_max], np.int16)
bin_centers = list()
bin_counts = list()

# 1d analysis

# binning just for theta_z
for i in range(max_bins):
    val = min_val + delta * (i + 0.5)
    # Determine which configurations lie in this bin.
    in_bin = (val-delta/2 <= th_ik[indices]) & (th_ik[indices] < val+delta/2) 
  
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

target_temp = 372.2
target_beta = 1/(kB*target_temp)
u_kn = np.zeros([K, N_max])
# populate diagonal blocks in MBAR array
for k in range(K):
    u_kn[k] = target_beta * uo_ik[k]
u_n = np.reshape(u_kn, N_sims*N_max)

[f_i, df_i] = my_mbar.computePMF(u_kn, bin_kn, nbins)
bin_centers = np.array(bin_centers)

prob_i = np.exp(-f_i)
theta_axis = bin_centers

plt.figure()
plt.plot(theta_axis, f_i, color='#006d2c', marker='o')
plt.fill_between(theta_axis, f_i - 2*df_i, f_i+2*df_i, color="#006d2c", alpha=0.4)
plt.xlim(-1.0, 0.0)
# plt.ylim(-0.2, 20)
plt.xlabel(r'$\theta$', fontsize=28)
plt.ylabel(r'$F(\theta)$', fontsize=28)
plt.show()

# plot E as a function of theta_z bins
th_flattened = th_ik.reshape(-1)
uo_flattened = uo_ik.reshape(-1)
E_bin, edges, y = sp.stats.binned_statistic(th_flattened, uo_flattened, statistic='median', bins=100)
# new_bins = 0.5 * (edges[1:] + edges[:-1])
# plt.plot(new_bins, E_bin, 'r', linewidth=5, alpha=0.4)
# plt.plot(new_bins, E_bin, 'ro', markersize=9, label='binned energies from all replicas')
# plt.xlabel(r'$\theta_z$', fontsize=28)
# plt.ylabel(r'$\langle E\rangle$', fontsize=28)
# plt.legend(loc='best', fontsize=28)
# plt.show()

# code to get an estimate of transition temeprature based on current dF and dE values
# dF already in units of beta*F

disord_min = np.mean(theta_axis[np.where(np.abs(theta_axis - -0.30) < 0.05)])
ord_min = np.mean(theta_axis[np.where(np.abs(theta_axis - -0.725) < 0.01)])
disord_F = np.mean(f_i[np.where(np.abs(theta_axis - disord_min) < 0.05)])
ord_F = np.mean(f_i[np.where(np.abs(theta_axis - ord_min) < 0.01)])

disord_E = np.mean( E_bin[ np.where(np.abs(new_bins - disord_min) < 0.05) ] )
ord_E = np.mean( E_bin[ np.where(np.abs(new_bins - ord_min) < 0.05) ] )

dF = ord_F - disord_F
dE = ord_E - disord_E
beta_transition = target_beta - dF / dE
temp_transition = 1 / (kB * beta_transition)
print "ord f = {} disord f = {} ord E = {} disord E = {} transition temp = {}K".format(ord_F, disord_F, ord_E, disord_E, temp_transition)

# # fit gaussians to free energy minima
# def gauss(x, x0, v, a):
#     return a * np.exp(-(x - x0)**2 / (2 * v))
# 
# ord_spring = np.where((theta_axis >= -0.75) * (theta_axis <=-0.69))
# x_ord = theta_axis[ord_spring]
# y_ord = prob_i[ord_spring]
# # popt, pcov = curve_fit(gauss, x_ord, y_ord, p0=(-0.75, 0.25, 0.8))
# popt, pcov = curve_fit(gauss, x_ord, y_ord)
# fit_ord = gauss(x_ord, popt[0], popt[1], popt[2])
# # print popt[0], popt[1], popt[2]
# 
# disord_spring = np.where((theta_axis >= -0.28) * (theta_axis <=0.15))
# x_disord = theta_axis[disord_spring]
# y_disord = prob_i[disord_spring]
# popt, pcov = curve_fit(gauss, x_disord, y_disord)
# fit_disord = gauss(x_disord, popt[0], popt[1], popt[2])
# # print popt[0], popt[1], popt[2]

# compute probability area of different phases and plot probability distribution
ord_indices = np.where(theta_axis <= -0.53)
disord_indices = np.where(theta_axis > -0.53)

area_ord = integrate.simps(prob_i[ord_indices], theta_axis[ord_indices])
area_disord = integrate.simps(prob_i[disord_indices], theta_axis[disord_indices])
print area_ord, area_disord

plt.figure()
plt.plot(theta_axis, prob_i, color='#006d2c', marker='o')
plt.xlabel(r'$\theta_z$', fontsize=28)
plt.ylabel(r'$P(\theta_z)$', fontsize=28)
plt.show()

