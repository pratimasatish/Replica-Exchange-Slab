import pymbar
import argparse
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from matplotlib import cm
from scipy.optimize import curve_fit
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D

# parser = argparse.ArgumentParser(description="")
# parser.add_argument("-dim", type=float, default=1, help="order parameter dimensionality [1d = th_z, 2d = th_z and th_x]")
# args = parser.parse_args()

# initialise list of temperatures 
kB = 1.3806503 * 6.0221415 / 4184.0

namelist_1 = np.arange(-0.8500, -0.0240, 0.0250)
namelist_2 = np.arange(0.0, 0.1040, 0.0250)
namelist = np.concatenate(( namelist_1, namelist_2 ))
N_bias = len(namelist)
full_namelist = np.concatenate(( namelist, namelist, namelist ))
T_bias = 5000

temp_list = np.concatenate(( np.ones(N_bias) * 350.18, np.ones(N_bias) * 355.0, np.ones(N_bias) * 350.0 ))
beta_list = 1/(kB * temp_list)

k_list = np.ones(N_bias * 3) * 15000.0

N_sims = len(temp_list)
theta_ik = []
theta_x_ik = []
VO_ik = []
UO_ik = []
 
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=28)
plt.rc('font', weight='bold')
plt.rc('xtick', labelsize=28)
plt.rc('ytick', labelsize=28)

# build up theta matrix

for biasval in full_namelist[:N_bias]:
    data = np.genfromtxt('theta-350.18.{:1.4f}.txt'.format(biasval))
    data = data.reshape((-1, 240))
    data_i = np.mean(data, axis=1)
    theta_ik.append( data_i[2000:] )

for biasval in full_namelist[N_bias:2*N_bias]:
    data = np.genfromtxt('theta-355.{:1.4f}.txt'.format(biasval))
    data = data.reshape((-1, 240))
    data_i = np.mean(data, axis=1)
    theta_ik.append( data_i[2000:] )

for biasval in full_namelist[2*N_bias:3*N_bias]:
    data = np.genfromtxt('theta-350.{:1.4f}.txt'.format(biasval))
    data = data.reshape((-1, 240))
    data_i = np.mean(data, axis=1)
    theta_ik.append( data_i[2000:] )

# build up theta_x matrix

for biasval in full_namelist[:N_bias]:
    data = np.genfromtxt('theta_x-350.18.{:1.4f}.txt'.format(biasval))
    data = data.reshape((-1, 240))
    data_i = np.mean(data, axis=1)
    theta_x_ik.append( data_i[2000:] )

for biasval in full_namelist[N_bias:2*N_bias]:
    data = np.genfromtxt('theta_x-355.{:1.4f}.txt'.format(biasval))
    data = data.reshape((-1, 240))
    data_i = np.mean(data, axis=1)
    theta_x_ik.append( data_i[2000:] )

for biasval in full_namelist[2*N_bias:3*N_bias]:
    data = np.genfromtxt('theta_x-350.{:1.4f}.txt'.format(biasval))
    data = data.reshape((-1, 240))
    data_i = np.mean(data, axis=1)
    theta_x_ik.append( data_i[2000:] )

# build up potential matrix

for k, th in enumerate(namelist[:N_bias]):
    lines = np.genfromtxt("pot-350.18.{:1.4f}".format(th))
    VO_ik.append( lines[2000:] )
    dtheta_i = np.array(theta_ik[k]) - th
    UO_ik.append( lines[2000:] - 0.5 * k_list[k] * np.square(dtheta_i) )

for k, th in enumerate(namelist[:N_bias]):
    lines = np.genfromtxt("pot-355.{:1.4f}".format(th))
    VO_ik.append( lines[2000:] )
    dtheta_i = np.array(theta_ik[k + N_bias]) - th
    UO_ik.append( lines[2000:] - 0.5 * k_list[k + N_bias] * np.square(dtheta_i) )

for k, th in enumerate(namelist[:N_bias]):
    lines = np.genfromtxt("pot-350.{:1.4f}".format(th))
    VO_ik.append( lines[2000:] )
    dtheta_i = np.array(theta_ik[k + 2*N_bias]) - th
    UO_ik.append( lines[2000:] - 0.5 * k_list[k + 2*N_bias] * np.square(dtheta_i) )

N_k = [ len(VO_i) for VO_i in VO_ik ]
N_k = np.array(N_k)
u_mbar = np.zeros((len(VO_ik), sum(N_k)))
K = u_mbar.shape[0]
N = u_mbar.shape[1]

# make numpy arrays from data
N_max = max(N_k)
th_ik = np.zeros([K, N_max])
uo_ik = np.zeros([K, N_max])
th_x_ik = np.zeros([K, N_max])
k = 0
for line1, line2, line3 in zip(theta_ik, UO_ik, theta_x_ik):
    th_ik[k,0:len(line1)] = np.array(line1)
    uo_ik[k,0:len(line2)] = np.array(line2)
    th_x_ik[k,0:len(line1)] = np.array(line3)
    k = k + 1

# go row by row to evaluate configuration energy at each temperature
for k in range(K):
    # populate off-diagonal blocks in MBAR array; go column by column, i.e. config by config
    for i in range(N_k[k]):
        dtheta = th_ik[k, i] - full_namelist
        print k, i
        u_mbar[ :, sum(N_k[:k]) + i ] = beta_list * ( uo_ik[k,i] + 0.5 * k_list * np.square(dtheta) )

my_mbar = pymbar.MBAR(u_mbar, N_k)

target_temp = 353.5
target_beta = 1/(kB*target_temp)
u_kn = np.zeros([K, N_max])
# populate diagonal blocks in MBAR array
for k in range(K):
    u_kn[k] = target_beta * uo_ik[k]
u_n = np.reshape(u_kn, N_sims*N_max)
theta_n = [val for row in theta_ik for val in row]
theta_n = np.array(theta_n)
theta_x_n = [val for row in theta_x_ik for val in row]
theta_x_n = np.array(theta_x_n)

mask_kn = np.zeros([K,N_max], dtype=np.bool)
for k in range(0,K):
   mask_kn[k,0:N_k[k]] = True
# Create a list from this mask.
indices = np.where(mask_kn)
max_bins = 250
nbins = 0
min_val = theta_n.min()
max_val = theta_n.max()
min_val_x = theta_x_n.min()
max_val_x = theta_x_n.max()
delta = (max_val - min_val) / float(max_bins)
delta_x = (max_val_x - min_val_x) / float(max_bins)
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

[f_i, df_i] = my_mbar.computePMF(u_kn, bin_kn, nbins)
bin_centers = np.array(bin_centers)

prob_i = np.exp(-f_i)
thz = bin_centers

def gauss(x, x0, v, a):
    return a * np.exp(-(x - x0)**2 / (2 * v))

ord_spring = np.where((thz >= -0.75) * (thz <=-0.69))
x_ord = thz[ord_spring]
y_ord = prob_i[ord_spring]
# popt, pcov = curve_fit(gauss, x_ord, y_ord, p0=(-0.75, 0.25, 0.8))
popt, pcov = curve_fit(gauss, x_ord, y_ord)
fit_ord = gauss(x_ord, popt[0], popt[1], popt[2])
print popt[0], popt[1], popt[2]

disord_spring = np.where((thz >= -0.28) * (thz <=0.15))
x_disord = thz[disord_spring]
y_disord = prob_i[disord_spring]
popt, pcov = curve_fit(gauss, x_disord, y_disord)
fit_disord = gauss(x_disord, popt[0], popt[1], popt[2])
print popt[0], popt[1], popt[2]

# compute and plot PMF as function of theta_z

plt.figure()
plt.fill_between(thz, f_i-2*df_i, f_i+2*df_i, color="#2020CC", alpha=0.4, linewidth=3)
plt.plot(thz, f_i, 'bo')
plt.xlabel(r'$\langle\theta_z\rangle$', fontsize=40)
plt.ylabel(r'$-\beta F$', fontsize=40)
plt.savefig('PMF-{}K.png'.format(target_temp))

ord_indices = np.where(thz <= -0.557)
disord_indices = np.where(thz > -0.557)

area_ord = integrate.simps(prob_i[ord_indices], thz[ord_indices])
area_disord = integrate.simps(prob_i[disord_indices], thz[disord_indices])
print area_ord, area_disord

# # 2d analysis
# 
# # binning for theta_z as well as theta_x
# for i in range(max_bins):
#     for j in range(max_bins):
#         val = min_val + delta * (i + 0.5)
#         val_x = min_val_x + delta_x * (j + 0.5)
#         # Determine which configurations lie in this bin.
#         in_bin = (val-delta/2 <= th_ik[indices]) & (th_ik[indices] < val+delta/2) & (val_x-delta_x/2 <= th_x_ik[indices]) & (th_x_ik[indices] < val_x+delta_x/2)
#       
#         # Count number of configurations in this bin.
#         bin_count = in_bin.sum()
#     
#         # Generate list of indices in bin.
#         indices_in_bin = (indices[0][in_bin], indices[1][in_bin])
#     
#         if (bin_count > 0):
#             bin_centers.append( (val, val_x) )
#             bin_counts.append( bin_count )
#      
#             # assign these conformations to the bin index
#             bin_kn[indices_in_bin] = nbins
#      
#             # increment number of bins
#             nbins += 1
# 
# [f_i, df_i] = my_mbar.computePMF(u_kn, bin_kn, nbins)
# bin_centers = np.array(bin_centers)
# 
# thz = bin_centers[:,0]
# thx = bin_centers[:,1]
# thzmin = thz.min()
# thzmax = thz.max()
# thxmin = thx.min()
# thxmax = thx.max()
# 
# thznew = np.linspace(thzmin + 0.01, thzmax - 0.01, 100)
# thxnew = np.linspace(thxmin + 0.01, thxmax - 0.03, 100)
# func = interpolate.bisplrep(thz, thx, f_i)
# fnew = interpolate.bisplev(thznew, thxnew, func)
# fnew = fnew.reshape(-1, 100)
# fnew = fnew.T
# Z, X = np.meshgrid(thznew, thxnew)
# 
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(Z, X, fnew, cmap=cm.hot, linewidth=0, antialiased=False, alpha=0.3)
# ax.scatter(thz, thx, f_i, c='k', alpha=1.0)
# plt.show()


