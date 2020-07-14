import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.mplot3d import Axes3D
import math

parser = argparse.ArgumentParser(description="")
parser.add_argument("-bias", type=str, help="bias value to analyse")
parser.add_argument("-steps", type=int, default=5000, help="number of time steps to analyse")
parser.add_argument("-pbc", action='store_true', help="whether to apply PBCs or not")
parser.add_argument("-savefigs", action='store_true', help="whether to save figures of CG lattice to make movie or not")
# parser.add_argument("-remove_NN", action='store_true', help="whether to remove sites with less than 4 NN's")
args = parser.parse_args()

lwidth = 4.0
plt.rc('text', usetex=True, fontsize=30)
rcParams['text.latex.preamble'] = [r'\usepackage{helvet} \usepackage{sfmath}', r'\usepackage{upgreek}' ]
rcParams['axes.linewidth'] = 1.0*lwidth
rcParams['xtick.major.width'] = 1.0*lwidth
rcParams['xtick.major.size']  = 2.0*lwidth
rcParams['ytick.major.width'] = 1.0*lwidth
rcParams['ytick.major.size']  = 2.0*lwidth
plt.rc('lines', linewidth=4)
plt.rc('legend', frameon=False)

data = np.genfromtxt('theta' + args.bias + '.txt', delimiter=' ')
data = data.reshape((-1,20,12))
Lx = 82.4293
Lz = 81.004

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
theta_t = theta_lat.reshape((-1, 240))
theta_lat = theta_lat.reshape((-1, 20, 12))
# print theta_lat.shape, theta_t.shape

# plt.clf();plt.hist(theta_lat[:,:].flatten(), bins=100, normed=True);plt.show()
# print theta_lat[:,:].flatten().var()

theta_mean = np.mean(theta_lat, axis=0)

# lattice plot of full data time averaged means

# print theta_mean.shape
# print 'x-mean', np.mean(theta_mean, axis=0)
# print 'z-mean', np.mean(theta_mean, axis=1)
plt.figure(figsize=(12,9))
plt.imshow(theta_mean.T, aspect=1/0.6, cmap="seismic_r", origin="lower", interpolation="none", vmin=-0.9, vmax=-0.6)
# plt.imshow(theta_mean, aspect=0.6, cmap="seismic_r", origin="lower", interpolation="none", vmin=-0.8, vmax=-0.1)
plt.yticks(np.arange(0, 12, 1))
plt.xticks(np.arange(0, 20, 1))
plt.ylim(-0.5,11.5)
plt.xlim(-0.5,19.5)
for i in np.arange(-0.5,12,1.0):
    plt.hlines(i, -0.5, 19.5, linestyle='solid', linewidth=2)
for i in np.arange(-0.5,19,1.0):
    plt.vlines(i, -0.5, 11.5, linestyle='solid', linewidth=2)
plt.xlabel('X')
plt.ylabel('Z')
plt.colorbar()
plt.savefig('lattice-{}.pdf'.format(args.bias))
plt.clf()

avg =  theta_lat.mean()
theta_mean = np.mean(theta_lat, axis=(1,2))
print "var = {} mean = {}".format(theta_mean.var(), theta_mean.mean())
# theta_mean -= avg
# print theta_mean.shape
# theta_lat -= avg 

if args.savefigs:
    name_arr = range(0, theta_lat.shape[0], 10)
    name_arr = np.array(name_arr)
    for j in range(len(name_arr)):
    # for j in range(0, 100, 20):
        matr = theta_lat[name_arr[j]].transpose()
        plt.imshow(matr, aspect=1.7, cmap="seismic_r", origin="lower", interpolation="none", vmin=-0.4, vmax=0.1)
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


# spatial correlation for full 100 ns
all_corr_xz = []
all_cov_xz = []
samplez = 10
# DT = t_steps / samplez
X = theta_lat.shape[1]
Z = theta_lat.shape[2]
DT = theta_lat.shape[0] / samplez
theta_lat -= np.mean(theta_lat)
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
    all_cov_xz.append(cov_xz)

all_corr_xz = np.array(all_corr_xz)
all_cov_xz = np.array(all_cov_xz)
m_corr_xz = np.mean(all_corr_xz, axis=0)
d_corr_xz = np.std (all_corr_xz, axis=0) / np.sqrt(samplez)

# correlation plots with first data point removed
plt.figure(figsize=(14,9))
plt.plot(range(0, Z/2), m_corr_xz[0,:], c='#01665e', label=r"$\textrm{Z}$", lw=4, marker='o', markersize=9)
plt.fill_between(range(0, Z/2), m_corr_xz[0,:] - d_corr_xz[0,:], m_corr_xz[0,:] + d_corr_xz[0,:], color='#01665e', alpha=0.4)
plt.plot(range(0, X/2), m_corr_xz[:,0], c='#8c510a', label=r"$\textrm{X}$", lw=4, marker='o', markersize=9)
plt.fill_between(range(0, X/2), m_corr_xz[:,0] - d_corr_xz[:,0], m_corr_xz[:,0] + d_corr_xz[:,0], color='#8c510a', alpha=0.4)
# plt.errorbar(range(1, X/2), m_corr_xz[1:,0], d_corr_xz[1:,0], c='b', label="X", linewidth=2)
# plt.errorbar(range(0, X/2), m_corr_xz[:,0], d_corr_xz[:,0], c='b', marker='o', label="X", linewidth=2)
plt.hlines(0, -0.1, X/2, linestyles="dashed")
plt.xlim(0.,X/2)
plt.ylim(-1.05,1.05)
# plt.errorbar(range(1, Z/2), m_corr_xz[0,1:], d_corr_xz[0,1:], c='g', label="Z", linewidth=2)
# plt.errorbar(range(0, Z/2), m_corr_xz[0,:], d_corr_xz[0,:], c='g', marker='v', label="Z", linewidth=2)
plt.xlabel(r'$\mathsf{r}$')
plt.ylabel(r'$\mathsf{G}_\theta\mathsf{(r)}$')
plt.legend(loc='upper right')
plt.savefig('solv-corr-xz-{}K.pdf'.format(args.bias))
plt.clf()

theta_lat = []
for i in range(T):
    theta_lat.append( ligdata[i, :, :].flatten() )
theta_lat = np.array(theta_lat)

# make coarse grained lattice
cd_data = np.genfromtxt('paired-sites.txt')
n_cd = len(cd_data)
# print n_cd

# blocks of four (average in x- and z-directions)
cg_lat = []

for i in range(n_cd):
    cg_lat.append((cd_data[i][2] + 4.12*0.5, cd_data[i][3], cd_data[i][4] + 6.75*0.5))

cg_lat = np.array(cg_lat)
nn_dist = np.sqrt( (4.12*0.5) **2 + (6.75*0.5)**2 )

# blocks of two (average only in x-direction)
cg_lat = []

for i in range(n_cd):
    cg_lat.append((cd_data[i][2] + 4.12*0.5, cd_data[i][3], cd_data[i][4]))

cg_lat = np.array(cg_lat)
nn_dist = 4.12*0.5

cg_xz = np.transpose( np.array((cg_lat[:, 0], cg_lat[:, 2])) )
# print cg_xz
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
print cg_theta.shape

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
        plt.savefig('lat-' + args.bias + '-{:05d}.png'.format(j))
#         plt.savefig('lat-0.7250-{:05d}.png'.format(j))
        plt.clf()

th_av = np.mean(cg_theta)
th_std = np.std(cg_theta)
# beta = 1/(400 * 1.3806503 * 6.0221415 / 4184.0)
# bins = np.linspace(-1.0, 0.5, 100)
# plt.hist(cg_theta[:,:,:].flatten(), bins=bins, normed=True, histtype='stepfilled', alpha=0.5)
# plt.plot(bins, np.exp(-beta * (bins - th_av)**2 / (2 * th_std * th_std) ) / np.sqrt(2 * np.pi * th_std *th_std), linewidth=3, color='b')
# plt.show()

theta_lat = theta_lat.reshape((-1, 20, 12))
print theta_lat.shape

# joint probability distributions of coarse-grained angle for a specific value of x
X = cg_mean.shape[0]
Z = cg_mean.shape[1]

plt.figure(figsize=(12,9))
plt.imshow(cg_mean.T, aspect=1/0.6, cmap="PRGn_r", origin="lower", interpolation="none", vmin=-0.9, vmax=-0.1)
plt.yticks(np.arange(0, 12, 1))
plt.xticks(np.arange(0, 20, 1))
plt.ylim(-0.5,11.5)
plt.xlim(-0.5,19.5)
for i in np.arange(-0.5,12,1.0):
    plt.hlines(i, -0.5, 19.5, linestyle='solid', linewidth=2)
for i in np.arange(-0.5,19,1.0):
    plt.vlines(i, -0.5, 11.5, linestyle='solid', linewidth=2)
plt.ylabel('Z', fontsize=28)
plt.xlabel('X', fontsize=28)
plt.colorbar()
plt.savefig('CG-lattice-{}.pdf'.format(args.bias))
plt.clf()

# # block averages for 25ns chunks
# 
# # first non-CG lattice
# 
# plt.clf()
# # plot time-avg of non-CG lattice over 25ns
# plt.figure(figsize=(9,9))
# plt.imshow(np.mean(theta_lat[0:1250], axis=0).T, aspect=1/0.6, cmap="seismic_r", origin="lower", interpolation="none", vmin=-0.15, vmax=0.00)
# plt.yticks(np.arange(0, 12, 1))
# plt.xticks(np.arange(0, 20, 1))
# plt.ylim(-0.5,11.5)
# plt.xlim(-0.5,19.5)
# for i in np.arange(-0.5,12,1.0):
#     plt.hlines(i, -0.5, 19.5, linestyle='solid', linewidth=2)
# for i in np.arange(-0.5,19,1.0):
#     plt.vlines(i, -0.5, 11.5, linestyle='solid', linewidth=2)
# plt.colorbar()
# plt.savefig('25-1.png')
# # plt.show()
# 
# plt.clf()
# # plot time-avg of non-CG lattice over second 25ns
# plt.figure(figsize=(9,9))
# plt.imshow(np.mean(theta_lat[1250:2500], axis=0).T, aspect=1/0.6, cmap="seismic_r", origin="lower", interpolation="none", vmin=-0.15, vmax=0.00)
# plt.yticks(np.arange(0, 12, 1))
# plt.xticks(np.arange(0, 20, 1))
# plt.ylim(-0.5,11.5)
# plt.xlim(-0.5,19.5)
# for i in np.arange(-0.5,12,1.0):
#     plt.hlines(i, -0.5, 19.5, linestyle='solid', linewidth=2)
# for i in np.arange(-0.5,19,1.0):
#     plt.vlines(i, -0.5, 11.5, linestyle='solid', linewidth=2)
# plt.colorbar()
# plt.savefig('25-2.png')
# # plt.show()
# 
# plt.clf()
# # plot time-avg of non-CG lattice over third 25ns
# plt.figure(figsize=(9,9))
# plt.imshow(np.mean(theta_lat[2500:3750], axis=0).T, aspect=1/0.6, cmap="seismic_r", origin="lower", interpolation="none", vmin=-0.15, vmax=0.00)
# plt.yticks(np.arange(0, 12, 1))
# plt.xticks(np.arange(0, 20, 1))
# plt.ylim(-0.5,11.5)
# plt.xlim(-0.5,19.5)
# for i in np.arange(-0.5,12,1.0):
#     plt.hlines(i, -0.5, 19.5, linestyle='solid', linewidth=2)
# for i in np.arange(-0.5,19,1.0):
#     plt.vlines(i, -0.5, 11.5, linestyle='solid', linewidth=2)
# plt.colorbar()
# plt.savefig('25-3.png')
# # plt.show()
# 
# plt.clf()
# # plot time-avg of non-CG lattice over last 25ns
# plt.figure(figsize=(9,9))
# plt.imshow(np.mean(theta_lat[3750:5000], axis=0).T, aspect=1/0.6, cmap="seismic_r", origin="lower", interpolation="none", vmin=-0.15, vmax=0.00)
# plt.yticks(np.arange(0, 12, 1))
# plt.xticks(np.arange(0, 20, 1))
# plt.ylim(-0.5,11.5)
# plt.xlim(-0.5,19.5)
# for i in np.arange(-0.5,12,1.0):
#     plt.hlines(i, -0.5, 19.5, linestyle='solid', linewidth=2)
# for i in np.arange(-0.5,19,1.0):
#     plt.vlines(i, -0.5, 11.5, linestyle='solid', linewidth=2)
# plt.colorbar()
# plt.savefig('25-4.png')
# # plt.show()
# 
# plt.clf()
# # then CG lattice
# 
# # plot time-avg of CG lattice over 25ns
# plt.figure(figsize=(9,9))
# plt.imshow(np.mean(cg_theta[0:1250], axis=0).T, aspect=1/0.6, cmap="seismic_r", origin="lower", interpolation="none", vmin=-0.15, vmax=0.00)
# plt.yticks(np.arange(0, 12, 1))
# plt.xticks(np.arange(0, 20, 1))
# plt.ylim(-0.5,11.5)
# plt.xlim(-0.5,19.5)
# for i in np.arange(-0.5,12,1.0):
#     plt.hlines(i, -0.5, 19.5, linestyle='solid', linewidth=2)
# for i in np.arange(-0.5,19,1.0):
#     plt.vlines(i, -0.5, 11.5, linestyle='solid', linewidth=2)
# plt.colorbar()
# plt.savefig('cg-25-1.png')
# # plt.show()
# 
# plt.clf()
# # plot time-avg of CG lattice over second 25ns
# plt.figure(figsize=(9,9))
# plt.imshow(np.mean(cg_theta[1250:2500], axis=0).T, aspect=1/0.6, cmap="seismic_r", origin="lower", interpolation="none", vmin=-0.15, vmax=0.00)
# plt.yticks(np.arange(0, 12, 1))
# plt.xticks(np.arange(0, 20, 1))
# plt.ylim(-0.5,11.5)
# plt.xlim(-0.5,19.5)
# for i in np.arange(-0.5,12,1.0):
#     plt.hlines(i, -0.5, 19.5, linestyle='solid', linewidth=2)
# for i in np.arange(-0.5,19,1.0):
#     plt.vlines(i, -0.5, 11.5, linestyle='solid', linewidth=2)
# plt.colorbar()
# plt.savefig('cg-25-2.png')
# # plt.show()
# 
# plt.clf()
# # plot time-avg of CG lattice over third 25ns
# plt.figure(figsize=(9,9))
# plt.imshow(np.mean(cg_theta[2500:3750], axis=0).T, aspect=1/0.6, cmap="seismic_r", origin="lower", interpolation="none", vmin=-0.15, vmax=0.00)
# plt.yticks(np.arange(0, 12, 1))
# plt.xticks(np.arange(0, 20, 1))
# plt.ylim(-0.5,11.5)
# plt.xlim(-0.5,19.5)
# for i in np.arange(-0.5,12,1.0):
#     plt.hlines(i, -0.5, 19.5, linestyle='solid', linewidth=2)
# for i in np.arange(-0.5,19,1.0):
#     plt.vlines(i, -0.5, 11.5, linestyle='solid', linewidth=2)
# plt.colorbar()
# plt.savefig('cg-25-3.png')
# # plt.show()
# 
# plt.clf()
# # plot time-avg of CG lattice over last 25ns
# plt.figure(figsize=(9,9))
# plt.imshow(np.mean(cg_theta[3750:5000], axis=0).T, aspect=1/0.6, cmap="seismic_r", origin="lower", interpolation="none", vmin=-0.15, vmax=0.00)
# plt.yticks(np.arange(0, 12, 1))
# plt.xticks(np.arange(0, 20, 1))
# plt.ylim(-0.5,11.5)
# plt.xlim(-0.5,19.5)
# for i in np.arange(-0.5,12,1.0):
#     plt.hlines(i, -0.5, 19.5, linestyle='solid', linewidth=2)
# for i in np.arange(-0.5,19,1.0):
#     plt.vlines(i, -0.5, 11.5, linestyle='solid', linewidth=2)
# plt.colorbar()
# plt.savefig('cg-25-4.png')
# # plt.show()

mean = np.mean(cg_theta[:,:,:], axis=(0,1,2))
# cg_theta = cg_theta[::10, :, :]
cg_theta[:, :, :] -= mean

cg_all_corr_xz = []
cg_all_cov_xz = []
samplez = 20
# DT = t_steps / samplez
DT = cg_theta.shape[0] / samplez
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
    
    cg_all_corr_xz.append(corr_xz)
    cg_all_cov_xz.append(cov_xz)
#     ax.plot_surface(xv, zv, corr_xz.T)

cg_all_corr_xz = np.array(cg_all_corr_xz)
m_corr_xz = np.mean(cg_all_corr_xz, axis=0)
d_corr_xz = np.std(cg_all_corr_xz, axis=0) / np.sqrt(samplez)

# correlation plot in both directions
plt.figure(figsize=(14,9))
plt.plot(range(0, Z/2), m_corr_xz[0,:], c='#01665e', label=r"$\mathsf{Z}$", lw=4, marker='o', markersize=9)
plt.fill_between(range(0, Z/2), m_corr_xz[0,:] - d_corr_xz[0,:], m_corr_xz[0,:] + d_corr_xz[0,:], color='#01665e', alpha=0.4)
plt.plot(range(0, X/2), m_corr_xz[:,0], c='#8c510a', label=r"$\mathsf{X}$", lw=4, marker='o', markersize=9)
plt.fill_between(range(0, X/2), m_corr_xz[:,0] - d_corr_xz[:,0], m_corr_xz[:,0] + d_corr_xz[:,0], color='#8c510a', alpha=0.4)
# plt.errorbar(range(1, X/2), m_corr_xz[1:,0], d_corr_xz[1:,0], c='b', label="X", linewidth=2)
# plt.errorbar(range(0, X/2), m_corr_xz[:,0], d_corr_xz[:,0], c='b', marker='o', label="X", linewidth=2)
plt.hlines(0, -0.1, X/2, linestyles="dashed")
plt.xlim(-0.1,X/2)
plt.ylim(-0.2, 1.05)
# plt.errorbar(range(1, Z/2), m_corr_xz[0,1:], d_corr_xz[0,1:], c='g', label="Z", linewidth=2)
# plt.errorbar(range(0, Z/2), m_corr_xz[0,:], d_corr_xz[0,:], c='g', marker='v', label="Z", linewidth=2)
plt.xlabel(r'$\mathsf{r}$')
plt.ylabel(r'$\mathsf{G}_\phi\mathsf{(r)}$')
plt.legend(loc='upper right')
plt.savefig('solv-cg-corr-xz-{}K.pdf'.format(args.bias))
plt.clf()

# correlation plots with first data point removed
plt.figure(figsize=(12,9))
# plt.errorbar(range(1, X/2), m_corr_xz[1:,0], d_corr_xz[1:,0], c='b', label="X", linewidth=2)
# plt.errorbar(range(0, X/2), m_corr_xz[:,0], d_corr_xz[:,0], c='#8c510a', label=r"$\textrm{X}$", lw=4)
# plt.errorbar(np.linspace(0,4.12*(X/2-1), X/2), m_corr_xz[:,0], d_corr_xz[:,0], c='b', label=r"$\textbf{X}$", linewidth=2)
# plt.errorbar(np.linspace(4.12, X/2*4.12, X/2-1), m_corr_xz[1:,0], d_corr_xz[1:,0], c='b', label="X", linewidth=2)
# plt.hlines(0, 0, X/2, linestyles="dashed")
# plt.xlim([0.9,Lx/2])
plt.ylim(-0.1,1.05)
plt.xlim(0,Z/2)
# plt.errorbar(range(1, Z/2), m_corr_xz[0,1:], d_corr_xz[0,1:], c='g', label="Z", linewidth=2)
plt.plot(range(0, Z/2), m_corr_xz[0,:], c='#01665e', label=r"$\textrm{Z}$", lw=4, marker='o', markersize=9)
plt.fill_between(range(0, Z/2), m_corr_xz[0,:] - d_corr_xz[0,:], m_corr_xz[0,:] + d_corr_xz[0,:], color='#01665e', alpha=0.4)
# plt.errorbar(np.linspace(0,6.75*(Z/2-1), Z/2), m_corr_xz[0,:], d_corr_xz[0,:], c='g', label=r"$\textbf{Z}$", linewidth=2)
# plt.errorbar(np.linspace(6.75, Z/2*6.75, Z/2-1), m_corr_xz[0,1:], d_corr_xz[0,1:], c='g', label="Z", linewidth=2)
plt.hlines(0, -0.1, Z/2, linestyles="dashed")
plt.xlabel(r'$\textrm{z}$')
plt.ylabel(r'$\textrm{G}_\phi\textrm{(z)}$')
plt.savefig('solv-cg-corr-thetaz-{}K.pdf'.format(args.bias))

for i in range(X/2):
    for j in range(Z/2):
        print "{} {} {} {}".format(i, j, m_corr_xz[i,j], d_corr_xz[i,j])

