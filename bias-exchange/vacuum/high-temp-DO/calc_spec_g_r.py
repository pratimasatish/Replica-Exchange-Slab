import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("-chain", type=int, default = 18, help="chain length")
parser.add_argument("-L", type=int, default = 81, help="box dimensions")
parser.add_argument("-N", type=int, default = 4320, help="# of particles")
args = parser.parse_args()

# data = np.genfromtxt('test-data.xyz')
data = np.genfromtxt('lig.last1250')
# RESHAPE TO TIME x PARTICLES
n_particle = args.N
# data_t = data.reshape((-1, n_particle))
T = int(len(data) / n_particle)
# print T
Lx = args.L
Lz = args.L
Ly = max(data[:,2]) - min(data[:,2])
rho = n_particle / (Lx * Ly * Lz)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=28)
plt.rc('font', weight='bold')
plt.rc('xtick', labelsize=28)
plt.rc('ytick', labelsize=28)

nbins = 50
n_chain = args.chain
n_lig = n_particle / n_chain
# L_max is box diagonal
L_max = np.sqrt(Lx**2 + Ly**2)
bins = np.linspace(0, L_max/2, nbins)
delta = bins[1] - bins[0]
hist_m = np.zeros(nbins)
hist_bt = np.zeros(nbins)

# FOR MIDDLE LIGAND ATOMS (NON BONDED ATOMS ONLY)
print "#middle atom indices"
for i in range(T):
    print "timestep {}".format(i)
    for j in range(8, n_particle, n_chain):
        # STORE START OF CHAIN INDEX SO INDEX+N_CHAIN CAN BE REMOVED FROM CALCULATIONS
        chain_start = j - 8
        # MODIFY RANGE FOR k TO EXCLUDE CHAIN
        for k in range(chain_start + n_chain, n_particle):
#             if ( (k-8) % n_chain == 0 or (k-9) % n_chain == 0 ):
            index1 = i * n_particle + j
            index2 = i * n_particle + k
#               print index1, index2
            del_x = (data[index1][1] - data[index2][1])
            del_x -= Lx * round(del_x / Lx)
            del_y = (data[index1][2] - data[index2][2])
            del_z = (data[index1][3] - data[index2][3])
            del_z -= Lz * round(del_z / Lz)

            dist = np.sqrt( del_x**2 + del_y**2 + del_z**2 )
            bin_index = np.where( (dist - delta/2 <= bins) & (bins < dist + delta/2) )
#             print bin_index
            hist_m[bin_index] += 2

    for j in range(9, n_particle, n_chain):
        # STORE START OF CHAIN INDEX SO INDEX+N_CHAIN CAN BE REMOVED FROM CALCULATIONS
        chain_start = j - 9
        # MODIFY RANGE FOR k TO EXCLUDE CHAIN
        for k in range(chain_start + n_chain, n_particle):
#             if ( (k-8) % n_chain == 0 or (k-9) % n_chain == 0 ):
            index1 = i * n_particle + j
            index2 = i * n_particle + k
#               print index1, index2
            del_x = (data[index1][1] - data[index2][1])
            del_x -= Lx * round(del_x / Lx)
            del_y = (data[index1][2] - data[index2][2])
            del_z = (data[index1][3] - data[index2][3])
            del_z -= Lz * round(del_z / Lz)

            dist = np.sqrt( del_x**2 + del_y**2 + del_z**2 )
            bin_index = np.where( (dist - delta/2 <= bins) & (bins < dist + delta/2) )
#             print bin_index
            hist_m[bin_index] += 2

# NORMALISATION FOR MIDDLE ATOMS
vol = np.ones(nbins)
for i in range(1, nbins):
    if (bins[i] <= Ly / 2):
        # SPHERICAL SHELL
        vol[i] = 4 * np.pi * delta * bins[i] ** 2
    else:
        # SPHERICAL SHELL WITH NO CAPS
        vol[i] = 2 * np.pi * delta * bins[i] * Ly

hist_m /= ( vol * rho * T * (n_particle-1-n_chain) )

# PRINT OUT BINS AND HISTOGRAM VALUES
for i in range(nbins):
    print "{} {}".format(bins[i], hist_m[i])


# FOR TOP AND BOTTOM LIGAND ATOMS (NON BONDED ATOMS ONLY)
print "#top and bottom atom indices"
for i in range(T):
    print "timestep {}".format(i)
    for j in range(0, n_particle, n_chain):
        # STORE START OF CHAIN INDEX SO INDEX+N_CHAIN CAN BE REMOVED FROM CALCULATIONS
        chain_start = j
        # MODIFY RANGE FOR k TO EXCLUDE CHAIN
        for k in range(chain_start + n_chain, n_particle):
#             if ( k % n_chain == 0 or (k+1) % n_chain == 0 ):
            index1 = i * n_particle + j
            index2 = i * n_particle + k
#              print index1, index2
            del_x = (data[index1][1] - data[index2][1])
            del_x -= Lx * round(del_x / Lx)
            del_y = (data[index1][2] - data[index2][2])
            del_z = (data[index1][3] - data[index2][3])
            del_z -= Lz * round(del_z / Lz)

            dist = np.sqrt( del_x**2 + del_y**2 + del_z**2 )
            bin_index = np.where( (dist - delta/2 <= bins) & (bins < dist + delta/2) )
#             print bin_index
            hist_bt[bin_index] += 2

    for j in range(17, n_particle, n_chain):
        # STORE START OF CHAIN INDEX SO INDEX+N_CHAIN CAN BE REMOVED FROM CALCULATIONS
        chain_start = j - 17
        # MODIFY RANGE FOR k TO EXCLUDE CHAIN
        for k in range(chain_start + n_chain, n_particle):
#             if ( k % n_chain == 0 or (k+1) % n_chain == 0 ):
            index1 = i * n_particle + j
            index2 = i * n_particle + k
#              print index1, index2
            del_x = (data[index1][1] - data[index2][1])
            del_x -= Lx * round(del_x / Lx)
            del_y = (data[index1][2] - data[index2][2])
            del_z = (data[index1][3] - data[index2][3])
            del_z -= Lz * round(del_z / Lz)

            dist = np.sqrt( del_x**2 + del_y**2 + del_z**2 )
            bin_index = np.where( (dist - delta/2 <= bins) & (bins < dist + delta/2) )
#             print bin_index
            hist_bt[bin_index] += 2

# NORMALISATION FOR TOP AND BOTTOM ATOMS
vol = np.ones(nbins)
for i in range(1, nbins):
    if (bins[i] <= Ly):
        # SPHERICAL SHELL
        vol[i] = 2 * np.pi * delta * bins[i] ** 2
    else:
        # SPHERICAL SHELL WITH NO CAPS
        vol[i] = 2 * np.pi * delta * bins[i] * Ly

hist_bt /= ( vol * rho * T * (n_particle-1-n_chain) )

# PRINT OUT BINS AND HISTOGRAM VALUES
for i in range(nbins):
    print "{} {}".format(bins[i], hist_bt[i])

exit(1)

# PLOTTING
plt.plot(bins, hist, 'ro')
plt.plot(bins, hist, 'r', linewidth=4, alpha=0.5)
plt.show()



