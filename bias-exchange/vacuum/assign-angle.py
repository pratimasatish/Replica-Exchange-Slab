import numpy as np
import argparse
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math

parser = argparse.ArgumentParser(description="")
parser.add_argument("-bias", type=str, help="bias value to analyse")
parser.add_argument("-steps", type=int, default=5000, help="number of time steps to analyse")
parser.add_argument("-pbc", action='store_true', help="whether to apply PBCs or not")
parser.add_argument("-savefigs", action='store_true', help="whether to save figures of CG lattice to make movie or not")
# parser.add_argument("-remove_NN", action='store_true', help="whether to remove sites with less than 4 NN's")
args = parser.parse_args()

data = np.genfromtxt('theta' + args.bias + '.txt', delimiter=' ')
data = data.reshape((-1,20,12))
Lx = 82.4293
Lz = 81.004

plt.rc('text', usetex=True)
plt.rc('font', family='serif', weight='bold')
plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
plt.rc('font', size=28)
plt.rc('xtick', labelsize=28)
plt.rc('ytick', labelsize=28)

data_txz = np.zeros((data.shape[0], data.shape[1], data.shape[2]))
data_txz[:, ::2, :] = data[:, 0:10, :]
data_txz[:, 1::2, :] = data[:, 10:20, :]
ligdata = data_txz
T = len(ligdata)

# get theta values at each lattice site as a function of time
theta_lat = []
for i in range(T):
    theta_lat.append( ligdata[i, :, :].flatten() )
theta_lat = np.array(theta_lat)
theta_lat = theta_lat.reshape((-1, 20, 12))
# print theta_lat.shape

# plt.clf();plt.hist(theta_lat[:,:].flatten(), bins=100, normed=True);plt.show()
# print theta_lat[:,:].flatten().var()

theta_mean = np.mean(theta_lat, axis=0)
# print theta_mean.shape
# print 'x-mean', np.mean(theta_mean, axis=0)
# print 'z-mean', np.mean(theta_mean, axis=1)
plt.imshow(theta_mean, aspect=0.6, cmap="seismic_r", origin="lower", interpolation="none", vmin=-0.4, vmax=-0.1)
# plt.imshow(theta_mean, aspect=0.6, cmap="seismic_r", origin="lower", interpolation="none", vmin=-0.8, vmax=-0.1)
plt.xticks(np.arange(0, 12, 1))
plt.yticks(np.arange(0, 20, 1))
plt.xlim(-0.5,11.5)
plt.ylim(-0.5,19.5)
for i in np.arange(-0.5,12,1.0):
    plt.vlines(i, -0.5, 19.5, linestyle='solid', linewidth=2)
for i in np.arange(-0.5,19,1.0):
    plt.hlines(i, -0.5, 11.5, linestyle='solid', linewidth=2)
plt.colorbar()
plt.show()

theta_lat -= theta_lat.mean()

if args.savefigs:
    name_arr = range(0, theta_lat.shape[0], 10)
    name_arr = np.array(name_arr)
    for j in range(len(name_arr)):
    # for j in range(0, 100, 20):
        matr = theta_lat[name_arr[j]].transpose()
        plt.imshow(matr, aspect=1.7, cmap="seismic_r", origin="lower", interpolation="none", vmin=-0.4, vmax=-0.1)
#         plt.imshow(matr, aspect=1.7, cmap="PRGn_r", origin="lower", interpolation="none", vmin=-0.8, vmax=-0.1)
        plt.yticks(np.arange(0, 12, 1))
        plt.xticks(np.arange(0, 20, 1))
        plt.ylim(-0.5,11.5)
        plt.xlim(-0.5,19.5)
        for i in np.arange(-0.5,12,1.0):
            plt.hlines(i, -0.5, 19.5, linestyle='solid', linewidth=2)
        for i in np.arange(-0.5,19,1.0):
            plt.vlines(i, -0.5, 11.5, linestyle='solid', linewidth=2)
        plt.colorbar()
        plt.savefig('lat-' + args.bias + '-{:05d}.png'.format(j))
#         plt.savefig('lat-0.7250-{:05d}.png'.format(j))
        plt.clf()


all_corr_xz = []
samplez = 10
# DT = t_steps / samplez
X = theta_lat.shape[1]
Z = theta_lat.shape[2]
DT = theta_lat.shape[0] / samplez
# fig = plt.figure()
# ax  = fig.add_subplot(111, projection='3d')
var_list = []
for sample in xrange(samplez):
    To = DT*sample
    Tf = DT*(sample+1)
    sub_txz = theta_lat[To:Tf, :, :]
    
    cov_xz = np.zeros((X/2, Z/2))
    for t in xrange(DT):
       for x in xrange(X/2):
            for z in xrange(Z/2):
#                 print cov_xz.shape
#                 print sub_txz[t,x,z].shape
#                 print sub_txz[t,x : x + X/2,z : z + Z/2].shape
                cov_xz += sub_txz[t, x, z] * sub_txz[t, x : x + X/2, z : z + Z/2]
    
    cov_xz /= (DT * X/2 * Z/2 )
    var_list.append(cov_xz[0,0])
    corr_xz = cov_xz / cov_xz[0,0]
#     corr_xz = cov_xz
    
    x = range(X/2)
    z = range(Z/2)
    xv, zv = np.meshgrid(x, z)
    
    all_corr_xz.append(corr_xz)

all_corr_xz = np.array(all_corr_xz)
m_corr_xz = np.mean(all_corr_xz, axis=0)
d_corr_xz = np.std (all_corr_xz, axis=0) / np.sqrt(samplez)
# print np.mean(np.array(var_list))
# print m_corr_xz[1:,0].shape, range(1,X/2)
 
# correlation plots with first data point removed
# plt.errorbar(range(1, X/2), m_corr_xz[1:,0], d_corr_xz[1:,0], c='b', label="X", linewidth=2)
plt.errorbar(range(0, X/2), m_corr_xz[:,0], d_corr_xz[:,0], c='b', label="X", linewidth=2)
plt.hlines(0, -1, X/2, linestyles="dashed")
plt.xlim([-1.0,X/2])
plt.ylim(-0.1,0.2)
# plt.errorbar(range(1, Z/2), m_corr_xz[0,1:], d_corr_xz[0,1:], c='g', label="Z", linewidth=2)
plt.errorbar(range(0, Z/2), m_corr_xz[0,:], d_corr_xz[0,:], c='g', label="Z", linewidth=2)
plt.hlines(0, -1, Z/2, linestyles="dashed")
plt.xlabel('x or z', fontsize=30)
plt.ylabel('G(x, z)', fontsize=30)
plt.legend(loc='upper right', fontsize=30)
plt.show()

# exit(1)

# alternate order parameter with differences instead of sums
T = theta_lat.shape[0]
X = theta_lat.shape[1]
Z = theta_lat.shape[2]

diff_angles = []

for i in range(T):
# only loop over even rows and take difference of odd row
    for x in range(0, X, 2):
        diff_angles.append(theta_lat[i][x] - theta_lat[i][x+1])
diff_angles = np.array(diff_angles)
# reshape for 10 row differences per time step to get a time-resolved array
diff_angles = diff_angles.reshape((-1, 10, 12))
diff_angles_time_mean = np.zeros(T)
for i in range(T):
    diff_angles_time_mean[i] = np.sum(diff_angles[i])/240

bins = np.linspace(diff_angles.min(), diff_angles.max(), 101)
bin_centres = 0.5 * (bins[1:] + bins[:-1])
pdist_diff, bins = np.histogram(diff_angles.flatten(), bins=bins, density=True)
pdist_diff_mean, bins = np.histogram(diff_angles_time_mean, bins=bins, density=True)

plt.plot(bin_centres, pdist_diff, 'ro', label='local data histogram')
plt.plot(bin_centres, pdist_diff, 'r', linewidth=4, alpha=0.7)
plt.plot(bin_centres, pdist_diff_mean, 'bv', label='global data histogram')
plt.plot(bin_centres, pdist_diff_mean, 'b', linewidth=4, alpha=0.7)
plt.legend(loc='best', fontsize=28)
plt.show()

# calculate correlation for this order parameter
all_corr_xz = []
samplez = 10
X = diff_angles.shape[1]
Z = diff_angles.shape[2]
DT = T/samplez
for sample in xrange(samplez):
    To = DT*sample
    Tf = DT*(sample+1)
    sub_txz = diff_angles[To:Tf, :, :]
    
    cov_xz = np.zeros((X/2, Z/2))
    for t in xrange(DT):
       for x in xrange(X/2):
            for z in xrange(Z/2):
                cov_xz += sub_txz[t, x, z] * sub_txz[t, x : x + X/2, z : z + Z/2]
 
    cov_xz /= (DT * X/2 * Z/2 )
    corr_xz = cov_xz / cov_xz[0,0]
    
    x = range(X/2)
    z = range(Z/2)
    xv, zv = np.meshgrid(x, z)
    
    all_corr_xz.append(corr_xz)
all_corr_xz = np.array(all_corr_xz)
m_corr_xz = np.mean(all_corr_xz, axis=0)
d_corr_xz = np.std (all_corr_xz, axis=0) / np.sqrt(samplez)

# plt.errorbar(range(1, X/2), m_corr_xz[1:,0], d_corr_xz[1:,0], c='b', label="X", linewidth=2)
plt.errorbar(range(0, X/2), m_corr_xz[:,0], d_corr_xz[:,0], c='b', label="X", linewidth=2)
# plt.errorbar(np.linspace(4.12, X/2*4.12, X/2-1), m_corr_xz[1:,0], d_corr_xz[1:,0], c='b', label="X", linewidth=2)
plt.hlines(0, 0, X/2, linestyles="dashed")
plt.xlim([0.0,X/2])
# plt.xlim([0.9,Lx/2])
plt.ylim(-0.1,1)
# plt.errorbar(range(1, Z/2), m_corr_xz[0,1:], d_corr_xz[0,1:], c='g', label="Z", linewidth=2)
plt.errorbar(range(0, Z/2), m_corr_xz[0,:], d_corr_xz[0,:], c='g', label="Z", linewidth=2)
# plt.errorbar(np.linspace(6.75, Z/2*6.75, Z/2-1), m_corr_xz[0,1:], d_corr_xz[0,1:], c='g', label="Z", linewidth=2)
plt.hlines(0, 0, Z/2, linestyles="dashed")
plt.xlabel('x or z', fontsize=30)
plt.ylabel('G(x, z)', fontsize=30)
plt.legend(loc='upper right', fontsize=30)
plt.show()

# exit(1)


theta_lat = []
for i in range(T):
    theta_lat.append( ligdata[i, :, :].flatten() )
theta_lat = np.array(theta_lat)


# make coarse grained lattice
cd_data = np.genfromtxt('paired-sites.txt')
n_cd = len(cd_data)
# print n_cd
cg_lat = []

for i in range(n_cd):
    cg_lat.append((cd_data[i][2] + 4.12*0.5, cd_data[i][3], cd_data[i][4] + 6.75*0.5))

cg_lat = np.array(cg_lat)
cg_xz = np.transpose( np.array((cg_lat[:, 0], cg_lat[:, 2])) )
# print cg_xz
nn_dist = np.sqrt( (4.12*0.5) **2 + (6.75*0.5)**2 )
# print cd_data[:,2].shape

pbc = True
full_indices = []
indices = []
# now loop over coarse grained lattice and get nearest neighbours for each point
for j in range(n_cd):
    x = cg_xz[j][0]
    z = cg_xz[j][1]
    # now find the Cd atoms closest to this coarse grained site
    dx = x - cd_data[:,2]
    dz = z - cd_data[:,4]
    # if args.pbc:
    if pbc:
        # apply PBCs
        dx = dx - Lx * np.round(dx / Lx)
        dz = dz - Lz * np.round(dz / Lz)

    dist = np.sqrt( dx**2 + dz**2 )
    full_indices.append( np.where( dist <= nn_dist+0.05 ) )
#     print j, indices

# remove_NN = args.remove_NN
remove_NN = None
if remove_NN == None:
    for j in range(n_cd):
        indices.append( full_indices[j] )
else:
    for j in range(n_cd):
        if ( np.array(full_indices[j]).shape[1] == 4 ):
            indices.append( full_indices[j] )        

# now assign theta values to each coarse grained site as a function of time
cg_theta = []
t_steps = theta_lat.shape[0]
# t_steps = args.steps
for i in range( T-t_steps, T ):
# for i in [T-1]:
    theta_t = []

    for j in range(len(indices)):
        theta_site = theta_lat[i-T + t_steps][indices[j]]
        theta_t.append(np.mean(theta_site))
    cg_theta.append(np.array(theta_t))

cg_theta = np.array(cg_theta)
if remove_NN == None:
    cg_theta = cg_theta.reshape((-1,20,12))
    cg_mean = np.mean(cg_theta, axis=0)
else:
    cg_mean = np.mean(cg_theta, axis=1)

# cg_mean = np.transpose(cg_mean)
# print cg_theta.shape

if args.savefigs:
    name_arr = range(0, t_steps, 10)
    name_arr = np.array(name_arr)
    for j in range(len(name_arr)):
    # for j in range(0, 100, 20):
        matr = cg_theta[name_arr[j]].transpose()
#         plt.imshow(matr, aspect=1.7, cmap="seismic_r", origin="lower", interpolation="none", vmin=-0.8, vmax=-0.1)
        plt.imshow(matr, aspect=1.7, cmap="PRGn_r", origin="lower", interpolation="none", vmin=-0.8, vmax=-0.1)
        plt.yticks(np.arange(0, 12, 1))
        plt.xticks(np.arange(0, 20, 1))
        plt.ylim(-0.5,11.5)
        plt.xlim(-0.5,19.5)
        for i in np.arange(-0.5,12,1.0):
            plt.hlines(i, -0.5, 19.5, linestyle='solid', linewidth=2)
        for i in np.arange(-0.5,19,1.0):
            plt.vlines(i, -0.5, 11.5, linestyle='solid', linewidth=2)
        plt.colorbar()
        plt.savefig('cg-lat-' + args.bias + '-{:05d}.png'.format(j))
#         plt.savefig('lat-0.7250-{:05d}.png'.format(j))
        plt.clf()

th_av = np.mean(cg_theta)
th_std = np.std(cg_theta)
beta = 1/(400 * 1.3806503 * 6.0221415 / 4184.0)
bins = np.linspace(-1.0, 0.5, 100)
plt.hist(cg_theta[:,:,:].flatten(), bins=bins, normed=True, histtype='stepfilled', alpha=0.7, color='green')
plt.xlabel(r'$\boldsymbol{\bar\theta}$', fontsize=28)
plt.ylabel(r'$\boldsymbol{P(\bar\theta)}$', fontsize=28)
plt.title(r'$\textbf{Probability distribution of 4-site average angle}$', fontsize=28)
# plt.plot(bins, np.exp(-beta * (bins - th_av)**2 / (2 * th_std * th_std) ) / np.sqrt(2 * np.pi * th_std *th_std), linewidth=3, color='b')
plt.show()

# print np.mean(np.mean(cg_mean, axis=1), axis=0)
# plt.imshow(cg_mean, aspect=0.6, cmap="seismic_r", origin="lower", interpolation="none", vmin=-0.8, vmax=-0.1)
plt.imshow(cg_mean, aspect=0.6, cmap="PRGn_r", origin="lower", interpolation="none", vmin=-0.8, vmax=-0.1)
plt.xticks(np.arange(0, 12, 1))
plt.yticks(np.arange(0, 20, 1))
plt.xlim(-0.5,11.5)
plt.ylim(-0.5,19.5)
for i in np.arange(-0.5,12,1.0):
    plt.vlines(i, -0.5, 19.5, linestyle='solid', linewidth=2)
for i in np.arange(-0.5,19,1.0):
    plt.hlines(i, -0.5, 11.5, linestyle='solid', linewidth=2)
plt.colorbar()
plt.show()

X = cg_mean.shape[0]
Z = cg_mean.shape[1]

mean = np.mean(cg_theta[:,:,:], axis=(0,1,2))
# print cg_theta.shape
# cg_theta = cg_theta[::10, :, :]
cg_theta[:, :, :] -= mean
# print cg_theta.shape

# print theta_lat[:,:].flatten().var()
# exit(1)

all_corr_xz = []
samplez = 20
# DT = t_steps / samplez
DT = cg_theta.shape[0] / samplez
fig = plt.figure()
ax  = fig.add_subplot(111, projection='3d')
dist_xz = np.zeros((X/2, Z/2))
for sample in xrange(samplez):
    To = DT*sample
    Tf = DT*(sample+1)
    sub_txz = cg_theta[To:Tf, :, :]
    
    cov_xz = np.zeros((X/2, Z/2))
    for t in xrange(DT):
       for x in xrange(X/2):
            for z in xrange(Z/2):
#                 print cov_xz.shape
#                 print sub_txz[t,x,z].shape
#                 print sub_txz[t,x : x + X/2,z : z + Z/2].shape
                cov_xz += sub_txz[t, x, z] * sub_txz[t, x : x + X/2, z : z + Z/2]
                dist_xz[x, z] = np.sqrt( (4.12*x)**2 + (6.75*z)**2 )
 
    cov_xz /= (DT * X/2 * Z/2 )
    corr_xz = cov_xz / cov_xz[0,0]
#     corr_xz = cov_xz
    
    x = range(X/2)
    z = range(Z/2)
    xv, zv = np.meshgrid(x, z)
    
    all_corr_xz.append(corr_xz)
#     ax.plot_surface(xv, zv, corr_xz.T)

# print dist_xz 
all_corr_xz = np.array(all_corr_xz)
m_corr_xz = np.mean(all_corr_xz, axis=0)
d_corr_xz = np.std (all_corr_xz, axis=0) / np.sqrt(samplez)
# print m_corr_xz[1:,0].shape, range(1,X/2)
ax.plot_surface(xv, zv, m_corr_xz.T)
plt.show()
for i in range(X/2):
    for j in range(Z/2):
        print i, j, m_corr_xz[i,j], d_corr_xz[i,j]

dist_r = np.unique(dist_xz)
corr_r = np.zeros(len(dist_r))
d_corr_r = np.zeros(len(dist_r))
for r in range(len(corr_r)):
    corr_r[r] =  np.mean(m_corr_xz[np.where(np.abs(dist_xz - dist_r[r]) <= 1e-4)])
#     d_corr_r[r] =  np.std(m_corr_xz[np.where(np.abs(dist_xz - dist_r[r]) <= 1e-4)])

plt.plot(dist_r, corr_r, 'b')
plt.plot(dist_r, corr_r, 'bo')
plt.xlabel('r', fontsize=30)
plt.ylabel('G(r)', fontsize=30)
plt.hlines(0, 1, Lz, linestyles="dashed")
plt.show()

# correlation plots with first data point removed
# plt.errorbar(range(1, X/2), m_corr_xz[1:,0], d_corr_xz[1:,0], c='b', label="X", linewidth=2)
# plt.errorbar(range(0, X/2), m_corr_xz[:,0], d_corr_xz[:,0], c='b', label=r"$\textbf{X}$", linewidth=2)
plt.errorbar(np.linspace(0,4.12*(X/2-1), X/2), m_corr_xz[:,0], d_corr_xz[:,0], c='b', label=r"$\textbf{X}$", linewidth=2)
# plt.errorbar(np.linspace(4.12, X/2*4.12, X/2-1), m_corr_xz[1:,0], d_corr_xz[1:,0], c='b', label="X", linewidth=2)
plt.hlines(0, 0, X/2*4.12, linestyles="dashed")
plt.xlim([0.0,X/2])
# plt.xlim([0.9,Lx/2])
plt.ylim(-0.1,1)
# plt.errorbar(range(1, Z/2), m_corr_xz[0,1:], d_corr_xz[0,1:], c='g', label="Z", linewidth=2)
# plt.errorbar(range(0, Z/2), m_corr_xz[0,:], d_corr_xz[0,:], c='g', label=r"$\textbf{Z}$", linewidth=2)
plt.errorbar(np.linspace(0,6.75*(Z/2-1), Z/2), m_corr_xz[0,:], d_corr_xz[0,:], c='g', label=r"$\textbf{Z}$", linewidth=2)
# plt.errorbar(np.linspace(6.75, Z/2*6.75, Z/2-1), m_corr_xz[0,1:], d_corr_xz[0,1:], c='g', label="Z", linewidth=2)
plt.hlines(0, 0, Z/2*6.75, linestyles="dashed")
plt.xlabel(r'\textbf{x or z}$', fontsize=30)
plt.ylabel(r'$\textbf{G(x, z)}$', fontsize=30)
plt.legend(loc='upper right', fontsize=30)
plt.show()


# pnt = ""
# for row in range(X):
#     strs = ["{:+1.4f}".format(x).zfill(5) for x in cg_mean[row, :]]
#     pnt = pnt + str(row).zfill(2) + " " + " ".join(strs)
#     pnt = pnt + "\n"
# print pnt

# print np.mean(cg_mean)
# print np.std(cg_mean)

# bins = np.linspace(-0.9, -0.5, 200)
# cg_hist, bins = np.histogram(cg_mean, bins=bins)
# cg_hist, bins = np.histogram(cg_theta[1000], bins=bins)
# bins = 0.5 * (bins[1:] + bins[:-1])
# plt.plot(bins, cg_hist, 'ro')
# plt.plot(bins, cg_hist)
# plt.show()


