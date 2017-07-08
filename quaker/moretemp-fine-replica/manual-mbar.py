import pymbar
import argparse
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from matplotlib import cm
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D

parser = argparse.ArgumentParser(description="")
parser.add_argument("-dim", type=float, default=1, help="order parameter dimensionality [1d = th_z, 2d = th_z and th_x]")
args = parser.parse_args()

# initialise list of temperatures 
kB = 1.3806503 * 6.0221415 / 4184.0
temp_list = [349.5, 350.0, 350.5, 351.0, 351.5, 352.0, 352.5, 353.0]
temp_list = np.array(temp_list)
beta_list = 1/(kB * temp_list)

N_sims = len(temp_list)
T = 2000
theta_ik = []
VO_ik = []
UO_ik = []
 
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

for k, temp in enumerate(temp_list):
    data = np.genfromtxt('theta{:3.1f}.txt'.format(temp))
    data = data.reshape((-1, 240))
    data_i = np.mean(data, axis=1)
    theta_ik.append( data_i[len(data_i) - T:] )

for k, temp in enumerate(temp_list):
    lines = np.genfromtxt("pot-new.{:3.1f}.txt".format(temp))
    VO_ik.append( lines[len(lines) - T:] )

N_k = [ len(VO_i) for VO_i in VO_ik ]
N_k = np.array(N_k)
u_mbar = np.zeros((len(VO_ik), sum(N_k)))
K = u_mbar.shape[0]
N = u_mbar.shape[1]

# make numpy arrays from data
N_max = max(N_k)
th_ik = np.zeros([K, N_max])
vo_ik = np.zeros([K, N_max])
k = 0
for line1, line2 in zip(theta_ik, VO_ik):
    th_ik[k,:] = np.array(line1)
    vo_ik[k,:] = np.array(line2)
    k = k + 1

# go row by row to evaluate configuration energy at each temperature
for k in range(K):
    # populate off-diagonal blocks in MBAR array; go column by column, i.e. config by config
    for i in range(N_k[k]):
        print k, i 
        u_mbar[ :, sum(N_k[:k]) + i ] = beta_list * ( vo_ik[k,i] )

my_mbar = pymbar.MBAR(u_mbar, N_k)

u_kn = []
target_temp = 350.18
target_beta = 1/(kB*target_temp)
# populate diagonal blocks in MBAR array
for i in range(K):
    u_kn.append(VO_ik[i] * target_beta)
u_kn = np.array(u_kn)
u_n = np.reshape(u_kn, N)

# one dimensional binning
nbins = 100
theta_n = [val for row in theta_ik for val in row]
theta_n = np.array(theta_n)
theta_n_sorted = np.sort(theta_n)
bins = np.append(theta_n_sorted[0::int(N/nbins)], theta_n_sorted.max()+0.005)
bin_widths = bins[1:] - bins[0:-1]
bin_n = np.zeros(theta_n.shape, np.int64)
bin_n = np.digitize(theta_n, bins) - 1

[f_i, df_i] = my_mbar.computePMF(u_n, bin_n, nbins)
f_i_corrected = f_i - np.log(bin_widths)
theta_axis = bins[:-1] * .5 + bins[1:] * .5

# two dimensional binning -- for some reason this binning works better probably because 
# the whole theta_z range is not uniformly sampled by replica exchange simulations
mask_kn = np.zeros([K,N_max], dtype=np.bool)
for k in range(0,K):
   mask_kn[k,0:N_k[k]] = True
# Create a list from this mask.
indices = np.where(mask_kn)
max_bins = 50
nbins = 0
min_val = -0.5
max_val = -0.1
delta = (max_val - min_val) / float(max_bins)
bin_kn = np.zeros([K,N_max], np.int16)
bin_centers = list()
bin_counts = list()
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

[f_i, df_i] = my_mbar.computePMF(u_kn, bin_kn, nbins)
theta_axis = np.array(bin_centers)

# compute and plot PMF as function of theta_z

prob_i = np.exp(-f_i)
plt.figure()
plt.plot(theta_axis, prob_i)
plt.show()

ord_indices = np.where(theta_axis < -0.557)
disord_indices = np.where(theta_axis > -0.557)

area_ord = integrate.simps(prob_i[ord_indices], theta_axis[ord_indices])
area_disord = integrate.simps(prob_i[disord_indices], theta_axis[disord_indices])
print area_ord, area_disord

