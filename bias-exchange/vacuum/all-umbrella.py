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

namelist_370int = np.array( [-0.750, -0.700, -0.690, -0.680, -0.670, -0.660, -0.650, -0.640, -0.630, -0.620, -0.610, -0.600, -0.575, -0.550] )
namelist_370mid = np.array( [-0.600, -0.580, -0.560, -0.540, -0.520, -0.500, -0.480, -0.460, -0.440, -0.420, -0.400, -0.380, -0.360, -0.340] )
namelist_370last = np.array( [-0.370, -0.350, -0.330, -0.310, -0.290, -0.270, -0.250, -0.230, -0.210, -0.190, -0.170, -0.150, -0.130, -0.110] )
namelist_380int = np.array( [-0.750, -0.700, -0.690, -0.680, -0.670, -0.660, -0.650, -0.640, -0.630, -0.620, -0.610, -0.600, -0.575, -0.550] )
namelist_380mid = np.array( [-0.600, -0.580, -0.560, -0.540, -0.520, -0.500, -0.480, -0.460, -0.440, -0.420, -0.400, -0.380, -0.360, -0.340] )
namelist_380last = np.array( [-0.370, -0.350, -0.330, -0.310, -0.290, -0.270, -0.250, -0.230, -0.210, -0.190, -0.170, -0.150, -0.130, -0.110] )

full_namelist = np.concatenate(( namelist_370int, namelist_370mid, namelist_370last, namelist_380int, namelist_380mid, namelist_380last ))
T_bias = 5000
N_bias = len(namelist_370int)

temp_list = np.concatenate(( np.ones(3 * N_bias) * 370.0, np.ones(3 * N_bias) * 380.0 ))
beta_list = 1/(kB * temp_list)
k_list = np.concatenate(( np.ones(N_bias) * 15000.0, np.ones(N_bias) * 7500.0, np.ones(N_bias) * 7500.0, np.ones(N_bias) * 7500.0, np.ones(N_bias) * 7500.0, np.ones(N_bias) * 7500.0 ))

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
for biasval in namelist_370int:
    data = np.genfromtxt('interface-370K/theta{:1.3f}.txt'.format(biasval))
    data = data.reshape((-1, 240))
    data_i = np.mean(data, axis=1)
    theta_ik.append( data_i )

for biasval in namelist_370mid:
    data = np.genfromtxt('interface-mid370K/theta{:1.3f}.txt'.format(biasval))
    data = data.reshape((-1, 240))
    data_i = np.mean(data, axis=1)
    theta_ik.append( data_i )

for biasval in namelist_370last:
    data = np.genfromtxt('interface-last370K/theta{:1.3f}.txt'.format(biasval))
    data = data.reshape((-1, 240))
    data_i = np.mean(data, axis=1)
    theta_ik.append( data_i )

for biasval in namelist_380int:
    data = np.genfromtxt('interface-start380K/theta{:1.3f}.txt'.format(biasval))
    data = data.reshape((-1, 240))
    data_i = np.mean(data, axis=1)
    theta_ik.append( data_i )

for biasval in namelist_380mid:
    data = np.genfromtxt('interface-mid380K/theta{:1.3f}.txt'.format(biasval))
    data = data.reshape((-1, 240))
    data_i = np.mean(data, axis=1)
    theta_ik.append( data_i )

for biasval in namelist_380last:
    data = np.genfromtxt('interface-last380K/theta{:1.3f}.txt'.format(biasval))
    data = data.reshape((-1, 240))
    data_i = np.mean(data, axis=1)
    theta_ik.append( data_i )

# build up potential matrix
for k, th in enumerate(namelist_370int):
    lines = np.genfromtxt("interface-370K/pot-new.{:1.3f}".format(th))
    VO_ik.append( lines[:] )
    dtheta_i = np.array(theta_ik[k]) - th
    UO_ik.append( lines[:] )

for k, th in enumerate(namelist_370mid):
    lines = np.genfromtxt("interface-mid370K/pot-new.{:1.3f}".format(th))
    VO_ik.append( lines[:] )
    dtheta_i = np.array(theta_ik[k + N_bias]) - th
    UO_ik.append( lines[:] )

for k, th in enumerate(namelist_370last):
    lines = np.genfromtxt("interface-last370K/pot-new.{:1.3f}".format(th))
    VO_ik.append( lines[:] )
    dtheta_i = np.array(theta_ik[k + 2*N_bias]) - th
    UO_ik.append( lines[:] )

for k, th in enumerate(namelist_380int):
    lines = np.genfromtxt("interface-start380K/pot-new.{:1.3f}".format(th))
    VO_ik.append( lines[:] )
    dtheta_i = np.array(theta_ik[k + 3*N_bias]) - th
    UO_ik.append( lines[:] )

for k, th in enumerate(namelist_380mid):
    lines = np.genfromtxt("interface-mid380K/pot-new.{:1.3f}".format(th))
    VO_ik.append( lines[:] )
    dtheta_i = np.array(theta_ik[k + 4*N_bias]) - th
    UO_ik.append( lines[:] )

for k, th in enumerate(namelist_380last):
    lines = np.genfromtxt("interface-last380K/pot-new.{:1.3f}".format(th))
    VO_ik.append( lines[:] )
    dtheta_i = np.array(theta_ik[k + 5*N_bias]) - th
    UO_ik.append( lines[:] )

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

target_temp = 375.6
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

disord_min = np.mean(theta_axis[np.where(np.abs(theta_axis - -0.30) < 0.05)])
ord_min = np.mean(theta_axis[np.where(np.abs(theta_axis - -0.725) < 0.01)])
disord_F = np.mean(f_i[np.where(np.abs(theta_axis - disord_min) < 0.05)])
ord_F = np.mean(f_i[np.where(np.abs(theta_axis - ord_min) < 0.01)])

dF = ord_F - disord_F
print "ord f = {} disord f = {}".format(ord_F, disord_F)

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

print "#beta*free energy at {}".format(target_temp)
for i in range(len(theta_axis)):
    print "{} {}".format(theta_axis[i], f_i[i])
