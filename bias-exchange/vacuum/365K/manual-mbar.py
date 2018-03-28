import pymbar
import argparse
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from matplotlib import cm
import scipy as sp
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D

# parser = argparse.ArgumentParser(description="")
# parser.add_argument("-dim", type=float, default=1, help="order parameter dimensionality [1d = th_z, 2d = th_z and th_x]")
# args = parser.parse_args()

kB = 1.3806503 * 6.0221415 / 4184.0
beta = 1/(kB * 365)

# namelist_1 = np.arange(0.25, 1.65, 0.15)
# namelist_2 = [1.50] 

# namelist = np.concatenate((namelist_1, namelist_2))
# namelist = np.arange(0.25, 1.65, 0.15)
namelist = [-0.750, -0.700, -0.675, -0.625, -0.600, -0.550, -0.500, -0.450, -0.425, -0.400, -0.375, -0.350, -0.275, -0.225, -0.175, -0.125]
namelist = np.array(namelist)

# initialise list of spring constants
k_list = np.ones(len(namelist)) * 1500.0

N_sims = len(namelist)
T = 5000
theta_ik = []
VO_ik = []
UO_ik = []
 
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('xtick', labelsize=28)
plt.rc('ytick', labelsize=28)

for k, biasval in enumerate(namelist):
#     data = np.genfromtxt('theta{:1.3f}.txt'.format(biasval))
    data = np.genfromtxt('theta{:1.3f}.txt'.format(biasval))
    data = data.reshape((-1, 240))
    data_t = np.mean(data, axis=1)
    theta_ik.append(data_t)
#      num_t = len(data_t) - 4131
#      theta_ik.append(data_t[num_t:])
#      theta_ik.append(data[::2])
# 
for k, th in enumerate(namelist):
    lines = np.genfromtxt("pot-new.{:1.3f}".format(th))
    VO_ik.append(lines[5000:])
    dtheta_i = np.array(theta_ik[k]) - th
    UO_ik.append( lines[5000:] - 0.5 * k_list[k] * np.square(dtheta_i) )
#     UO_ik.append( 0.5 * k_list[k] * np.square(dtheta_i) )

N_k = [ len(UO_i) for UO_i in UO_ik ]
N_k = np.array(N_k)
u_mbar = np.zeros((len(UO_ik), sum(N_k)))
K = u_mbar.shape[0]
N = u_mbar.shape[1]

ave_u = np.zeros(K)
for k in range(K):
    ave_u[k] = np.mean(UO_ik[k])
    print namelist[k], ave_u[k]

# exit(0)

# make numpy arrays from data
N_max = max(N_k)
th_ik = np.zeros([K, N_max])
uo_ik = np.zeros([K, N_max])
k = 0
for line1, line2 in zip(theta_ik, UO_ik):
    th_ik[k,:] = np.array(line1)
    uo_ik[k,:] = np.array(line2)
    k = k + 1

# th_ik = th_ik[:,::2]

# go row by row to evaluate configuration energy in each umbrella
for k in range(K):
    # populate off-diagonal blocks in MBAR array; go column by column, i.e. config by config
    for i in range(N_k[k]):
        dtheta = th_ik[k, i] - namelist		# deviation of current configuration from each umbrella centre 
#         print k, i 
#         u_mbar[ k, count:count+N_k[i] ] = beta * ( UO_ik[k] + 0.5 * k_list[k] * np.square(dtheta_k) ) - beta * ( UO_ik[i] + 0.5 * k_list[i] * np.square(dtheta_i) )
        u_mbar[ :, sum(N_k[:k]) + i ] = beta * ( uo_ik[k,i] + 0.5 * k_list * np.square(dtheta) )

my_mbar = pymbar.MBAR(u_mbar, N_k)

u_kn = []
# populate diagonal blocks in MBAR array
for i in range(K):
    u_kn.append(UO_ik[i] * beta)
u_kn = np.array(u_kn)
u_n = np.reshape(u_kn, N)
theta_n = [val for row in theta_ik for val in row]
theta_n = np.array(theta_n)

# one dimensional binning
# nbins = 250
# theta_n_sorted = np.sort(theta_n)
# bins = np.append(theta_n_sorted[0::int(N/nbins)], theta_n_sorted.max()+0.005)
# bin_widths = bins[1:] - bins[0:-1]
# bin_n = np.zeros(theta_n.shape, np.int64)
# bin_n = np.digitize(theta_n, bins) - 1
# 
# [f_i, df_i] = my_mbar.computePMF(u_n, bin_n, nbins)
# f_i_corrected = f_i - np.log(bin_widths)
# theta_axis = bins[:-1] * .5 + bins[1:] * .5

# two dimensional binning
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

# compute and plot PMF as function of theta_z

[f_i, df_i] = my_mbar.computePMF(u_kn, bin_kn, nbins)
theta_axis = np.array(bin_centers)

prob_i = np.exp(-f_i)
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
new_bins = 0.5 * (edges[1:] + edges[:-1])
# print 'E_bin shape = {}'.format(E_bin.shape)
# print 'theta_axis shape = {}'.format(theta_axis.shape)
plt.plot(new_bins, E_bin, 'r', linewidth=5, alpha=0.4)
plt.plot(new_bins, E_bin, 'ro', markersize=9, label='binned energies from all replicas')
plt.plot(namelist, ave_u, 'bv', markersize=10, label='averaged for each replica')
plt.xlabel(r'$\theta_z$', fontsize=28)
plt.ylabel(r'$\langle E\rangle$', fontsize=28)
plt.legend(loc='best', fontsize=28)
plt.show()

# print "free energies\n" 
# for j in range(nbins):
#     print theta_axis[j], f_i[j]
# 
# code to get an estimate of transition temeprature based on current dF and dE values
# dF already in units of beta*F

disord_min = np.mean(theta_axis[np.where(np.abs(theta_axis - -0.350) < 0.005)])
ord_min = np.mean(theta_axis[np.where(np.abs(theta_axis - -0.75) < 0.005)])
disord_F = np.mean(f_i[np.where(np.abs(theta_axis - disord_min) < 0.005)])
ord_F = np.mean(f_i[np.where(np.abs(theta_axis - ord_min) < 0.005)])

disord_E = np.mean( E_bin[ np.where(np.abs(new_bins - disord_min) < 0.01) ] )
ord_E = np.mean( E_bin[ np.where(np.abs(new_bins - ord_min) < 0.01) ] )

dF = ord_F - disord_F
dE = ord_E - disord_E
beta_transition = beta - dF / dE
temp_transition = 1 / (kB * beta_transition)
print "ord f = {} disord f = {} ord E = {} disord E = {} transition temp = {}K".format(ord_F, disord_F, ord_E, disord_E, temp_transition)

# code to subtract "straight line" from free energy profile
# if only chemical potential difference matters, assumption should be correct
# and subtraction should lead to coexistence (no chem pot difference!)
slope = (disord_F - ord_F) / (disord_min - ord_min)
line = slope * (theta_axis - ord_min) + ord_F
plt.plot(theta_axis, f_i - line, 'bs')
plt.show()

# get actual centres after applying bias
# can invert to get where to place simulation 
# centers to get preferred actual centres

th_cen = np.zeros(K)
for index in range(K):
    tot_en = f_i + 0.5 * k_list[index] * (theta_axis - namelist[index])**2
    abs_der = abs(np.gradient(tot_en, delta))
    th_arr = theta_axis[np.where( (abs_der - np.min(abs_der) ) <= 1.0e1)]
    th_cen[index] = np.mean(th_arr)
#     print namelist[index], th_cen[index]

plt.figure()
plt.plot(namelist, th_cen, 'bo')
plt.plot(namelist, th_cen, 'b', linewidth=4, alpha=0.5)
plt.plot(namelist, namelist, 'r', linewidth=4, alpha=0.5)
plt.plot(namelist, namelist, 'ro')
plt.xlabel(r'$\theta_0$', fontsize=28)
plt.ylabel(r'$\theta_c$', fontsize=28)
plt.show()

# plot all the biased histograms
new_uo = np.zeros((K, nbins))
for k in range(K):
    new_uo[k, :] = 0.5 * k_list[k] * (bin_centers - namelist[k])**2

color_list=iter(plt.cm.rainbow(np.linspace(0,1,K)))
for index in range(K):
    c = next(color_list)
    plt.hist(th_ik[index, :], bins=bin_centers, normed=True, label=r'$\theta = {}$'.format(namelist[index]), histtype='stepfilled', alpha=0.5, color=c)
#     plt.hist(th_ik[index, :], bins=bin_centers, label=r'$\theta = {}$'.format(namelist[index]), histtype='stepfilled', alpha=0.5)
    plt.plot(bin_centers, (np.sqrt(beta * k_list[k])) * np.exp(-beta * new_uo[index, :]), linewidth=3, color=c)
plt.xlabel(r'$\theta$', fontsize=28)
plt.ylabel(r'$P(\theta)$', fontsize=28)
plt.legend(loc='best', fontsize=20)
plt.show()

# reweight free energy to new target temperature
# target_temp = 372
# target_beta = 1/(kB * target_temp)
u_kn = []
# populate diagonal blocks in MBAR array
for i in range(K):
    u_kn.append(UO_ik[i] * beta_transition)
u_kn = np.array(u_kn)
u_n = np.reshape(u_kn, N)

[f_new, df_new] = my_mbar.computePMF(u_kn, bin_kn, nbins)

# also calculate a simple estimate of free energy at a new temperature
# compare to MBAR results
shifted_u = E_bin - np.mean(E_bin)
simple_f_new = f_i + (beta_transition - beta) * shifted_u
simple_f_new = simple_f_new - np.min(simple_f_new)

for i in range(len(new_bins)):
    print '{} {} {}'.format(new_bins[i], f_new[i], simple_f_new[i])

prob_new = np.exp(-f_new)
simple_prob_new = np.exp(-simple_f_new)
plt.figure()
plt.plot(theta_axis, f_new, color='#006d2c', marker='o', label='MBAR predicted free energy')
plt.fill_between(theta_axis, f_new - 2*df_new, f_new+2*df_new, color="#006d2c", alpha=0.4)
plt.plot(new_bins, simple_f_new, 'rv', markersize=8, label='simple free energy estimate')
plt.xlim(-1.0, 0.0)
plt.ylim(-2, 25)
plt.xlabel(r'$\theta$', fontsize=28)
plt.ylabel(r'$F_{\textrm{target}}(\theta)$', fontsize=28)
plt.title('free energy at {}$'.format(temp_transition), fontsize=32)
plt.legend(loc='best', fontsize=22)
plt.show()

# compute probability area of different phases at target temperature = original temperature
ord_indices = np.where(theta_axis < -0.575)
disord_indices = np.where(theta_axis > -0.575)

area_ord = integrate.simps(prob_new[ord_indices], theta_axis[ord_indices])
area_disord = integrate.simps(prob_new[disord_indices], theta_axis[disord_indices])

simple_area_ord = integrate.simps(simple_prob_new[ord_indices], theta_axis[ord_indices])
simple_area_disord = integrate.simps(simple_prob_new[disord_indices], theta_axis[disord_indices])

print area_ord, area_disord
print simple_area_ord, simple_area_disord

plt.figure()
plt.plot(theta_axis, prob_new, color='#006d2c', marker='o', label='MBAR predicted free energy')
plt.plot(new_bins, simple_prob_new, 'rv', label='simple free energy estimate')
plt.xlabel(r'$\theta_z$', fontsize=28)
plt.ylabel(r'$P(\theta_z)$', fontsize=28)
plt.legend(loc='best', fontsize=22)
plt.show()


