import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("-chain", type=int, default = 18, help="chain length")
parser.add_argument("-L", type=int, default = 81, help="box dimensions")
parser.add_argument("-N", type=int, default = 4320, help="# of particles")
args = parser.parse_args()

# data = np.genfromtxt('test-data.xyz')
data = np.genfromtxt('lig.1000last')
# RESHAPE TO TIME x PARTICLES
n_particle = args.N
# data_t = data.reshape((-1, n_particle))
T = int(len(data) / n_particle)
# print T
L_box = args.L
rho = float(n_particle) / (L_box)**3

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=28)
plt.rc('font', weight='bold')
plt.rc('xtick', labelsize=28)
plt.rc('ytick', labelsize=28)

nbins = 50
n_chain = args.chain
n_lig = n_particle / n_chain
bins = np.linspace(0, L_box/2, nbins)
delta = bins[1] - bins[0]
hist = np.zeros(nbins)

# FOR LIGAND MOLECULES (NON BONDED ATOMS ONLY)
for i in range(T):
#     print "timestep {}".format(i)
    for j in range(n_particle):
        # IF THIS IS START OF A CHAIN, STORE INDEX SO INDEX+N_CHAIN CAN BE REMOVED FROM CALCULATIONS
        if ( j % n_chain == 0):
            chain_start = j
        # MODIFY RANGE FOR k TO EXCLUDE CHAIN
        for k in range(chain_start + n_chain, n_particle):
            index1 = i * n_particle + j
            index2 = i * n_particle + k
#             print index1, index2
            del_x = (data[index1][1] - data[index2][1])
            del_x -= L_box * round(del_x / L_box)
            del_y = (data[index1][2] - data[index2][2])
            del_y -= L_box * round(del_y / L_box)
            del_z = (data[index1][3] - data[index2][3])
            del_z -= L_box * round(del_z / L_box)

            dist = np.sqrt( del_x**2 + del_y**2 + del_z**2 )
            bin_index = np.where( (dist - delta/2 <= bins) & (bins < dist + delta/2) )
#             print bin_index
            hist[bin_index] += 2


# FOR RANDOM PARTICLE IN A BOX (UNIFORM DISTRIBUTION)
# for i in range(T):
#     for j in range(n_particle):
#         for k in range(j+1, n_particle):
#             index1 = i * n_particle + j
#             index2 = i * n_particle + k
# #             print index1, index2
#             del_x = (data[index1][1] - data[index2][1])
#             del_x -= L_box * round(del_x / L_box)
#             del_y = (data[index1][2] - data[index2][2])
#             del_y -= L_box * round(del_y / L_box)
#             del_z = (data[index1][3] - data[index2][3])
#             del_z -= L_box * round(del_z / L_box)
# 
#             dist = np.sqrt( del_x**2 + del_y**2 + del_z**2 )
#             bin_index = np.where( (dist - delta/2 <= bins) & (bins < dist + delta/2) )
# #             print bin_index
#             hist[bin_index] += 2

# NORMALISATION ???????
vol = np.zeros(nbins)
vol[0] = (4 / 3) * np.pi * (delta/2) ** 3
# vol[1:] = (4 / 3) * np.pi * ( (bins[1:] + delta/2)**3 - (bins[1:] - delta/2)**3 )
vol[1:] = 4 * np.pi * bins[1:]**2 * delta
# hist /= ( vol * rho * T * (n_particle-1) )
hist /= ( vol * rho * T * (n_particle-1-n_chain) )

# PRINT OUT BINS AND HISTOGRAM VALUES
for i in range(nbins):
    print "{} {}".format(bins[i], hist[i])

# PLOTTING
plt.plot(bins, hist, 'ro')
plt.plot(bins, hist, 'r', linewidth=4, alpha=0.5)
plt.show()



