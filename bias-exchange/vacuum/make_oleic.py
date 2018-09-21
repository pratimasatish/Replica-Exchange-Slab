import numpy as np
import argparse
from random import randint

parser = argparse.ArgumentParser()
parser.add_argument("-nx", type=int, default=20, help="number of ligands in x-direction")
parser.add_argument("-nz", type=int, default=12, help="number of ligands in z-direction")
parser.add_argument("-cov", type=float, default=0.75, help="fraction of surface passivated by ligands")
parser.add_argument("-solv", action='store_true', help="solvate the system or not")
parser.add_argument("-fmt", type=str, choices=["xyz", "lmp"], default='xyz', help="output format")
args = parser.parse_args()

# data = np.genfromtxt('backbone.txt', delimiter=' ')
data = np.genfromtxt('from200-chain.txt', delimiter=' ')
# liq_hex = np.genfromtxt('liq_solv.txt', delimiter=' ')
liq_hex = np.genfromtxt('hex_solv.xyz', delimiter=' ')
x0 = data[0,1]
y0 = data[0,2]
z0 = data[0,3]
sx0 = np.min(liq_hex[:,1])
sy0 = np.min(liq_hex[:,2])
sz0 = np.min(liq_hex[:,3])
centred_data = np.zeros( (len(data), 3) )
centred_bot = np.zeros( (len(data), 3) )
centred_hex = np.zeros( (len(liq_hex), 3) )
nx = args.nx
nz = args.nz
# passivate at given% on both surfaces, periodic or randomly
n_molecules = nx * nz
n_top = int(args.cov * nx * nz)
n_bottom = int(0.5 * nx * nz)
n_chain = 18
n_bonds = n_chain - 1
n_angles = n_chain - 2
n_dihedrals = n_chain - 4

if args.solv:
    n_hexatoms = len(liq_hex) 
else:
    n_hexatoms = 0
n_hex = 6
n_hexmolecules = n_hexatoms / n_hex
n_hexbonds = n_hex - 1
n_hexangles = n_hex - 2
n_hexdihedrals = n_hex - 3
# xlo =  -5.0
# xhi =  80.0
# ylo = -60.0
# yhi = 120.0
# zlo =  -5.0
# zhi =  80.0

xlo =  0.0
xhi =  82.4293
ylo = -60.0
yhi = 120.0
zlo =  0.0
zhi =  81.004

rotate = np.array([ [np.cos(0.5 * np.pi), -np.sin(0.5 * np.pi), 0], [np.sin(0.5 * np.pi), np.cos(0.5 * np.pi), 0], [0, 0, 1.0] ])
reflect = np.array([ [1.0, 0, 0], [0, -1.0, 0], [0, 0, 1.0] ])

# centre the original backbone
for i in range(len(data)):
    centred_data[i,0] = data[i,1] - x0
    centred_data[i,1] = data[i,2] - y0
    centred_data[i,2] = data[i,3] - z0

# centre the hex chain
for i in range(len(liq_hex)):
    centred_hex[i,0] = liq_hex[i,1]
    centred_hex[i,1] = liq_hex[i,2]
    centred_hex[i,2] = liq_hex[i,3]

# for i in range(len(centred_data)):
#     centred_data[i, :] = rotate.dot(centred_data[i, :])

# reflected for the bottom surface
for i in range(len(centred_data)):
    centred_bot[i, :] = reflect.dot(centred_data[i, :])

if args.fmt == "xyz":
    print "{}".format(n_top * n_chain + n_bottom * n_chain + n_molecules * 16 + n_hexatoms)
    print "Comment Line"
    nmol = 0
    # atomic coordinates for top ligands here (given% random passivation)
    N_lig = int(args.cov * nx * nz)
    occ_indices = []

    while (N_lig > 0):
        site_x = randint(0, nx - 1)
        site_z = randint(0, nz - 1)
        site_num = [site_x, site_z]
        if site_num not in occ_indices:
            # place ligand
            coords = np.copy(centred_data)
            # shift molecule
            coords[:, 0] = coords[:, 0] + 4.12 * site_x
            coords[:, 2] = coords[:, 2] + 6.75 * site_z
#             print site_x, site_z        

            for k in range(n_chain):
                if (k == n_chain-1):
                    print "4       {:4.5f} {:4.5f} {:4.5f}".format( coords[k,0], coords[k,1], coords[k,2] )
                elif ((k == 7) or (k == 10)):
                    print "2       {:4.5f} {:4.5f} {:4.5f}".format( coords[k,0], coords[k,1], coords[k,2] )
                elif ((k == 8) or (k == 9)):
                    print "3       {:4.5f} {:4.5f} {:4.5f}".format( coords[k,0], coords[k,1], coords[k,2] )
                else:
                    print "1       {:4.5f} {:4.5f} {:4.5f}".format( coords[k,0], coords[k,1], coords[k,2] )
    
            nmol = nmol + 1
            occ_indices.append(site_num)
            N_lig = N_lig - 1

#     for i in range(nx):
#         for j in range(nz):
#             coords = np.copy(centred_data)
#             # shift molecule
#             coords[:, 0] = coords[:, 0] + 4.12 * i
#             coords[:, 2] = coords[:, 2] + 6.75 * j
#         
#             for k in range(n_chain):
#                 if (k == n_chain-1):
#                     print "4       {:4.5f} {:4.5f} {:4.5f}".format( coords[k,0], coords[k,1], coords[k,2] )
#                 elif ((k == 7) or (k == 10)):
#                     print "2       {:4.5f} {:4.5f} {:4.5f}".format( coords[k,0], coords[k,1], coords[k,2] )
#                 elif ((k == 8) or (k == 9)):
#                     print "3       {:4.5f} {:4.5f} {:4.5f}".format( coords[k,0], coords[k,1], coords[k,2] )
#                 else:
#                     print "1       {:4.5f} {:4.5f} {:4.5f}".format( coords[k,0], coords[k,1], coords[k,2] )
#     
#             nmol = nmol + 1

    # Cd-S surface coordinates
    for i in range(nx):
        for j in range(nz):
            coords = centred_data[0, :] - np.array([0.0, 1.56, 0.0])
            coords[0] = coords[0] + 4.12 * i
            coords[2] = coords[2] + 6.75 * j
            print "5	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[0], coords[1], coords[2] )
            print "6	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[0], coords[1], coords[2] + 2.5299 )
            print "5	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[0] - 2.06073, coords[1] - 1.18977, coords[2] + 3.375 )
            print "6	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[0] - 2.06073, coords[1] - 1.18977, coords[2] - 3.375 + 2.5299 )
            print "5	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[0] - 2.06073, coords[1] - 3.5693, coords[2] )
            print "6	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[0] - 2.06073, coords[1] - 3.5693, coords[2] + 2.5299 )
            print "5	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[0], coords[1] - 3.5693 - 1.18977, coords[2] + 3.375 )
            print "6	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[0], coords[1] - 3.5693 - 1.18977, coords[2] - 3.375 + 2.5299 )

            print "5	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[0], coords[1] - 7.13848, coords[2] )
            print "6	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[0], coords[1] - 7.13848, coords[2] + 2.5299 )
            print "5	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[0] - 2.06073, coords[1] - 7.13848 - 1.18977, coords[2] + 3.375 )
            print "6	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[0] - 2.06073, coords[1] - 7.13848 - 1.18977, coords[2] - 3.375 + 2.5299 )
            print "5	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[0] - 2.06073, coords[1] - 7.13848 - 3.5693, coords[2] )
            print "6	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[0] - 2.06073, coords[1] - 7.13848 - 3.5693, coords[2] + 2.5299 )
            print "5	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[0], coords[1] - 7.13848 - 3.5693 - 1.18977, coords[2] + 3.375 )
            print "6	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[0], coords[1] - 7.13848 - 3.5693 - 1.18977, coords[2] - 3.375 + 2.5299 )

    # bottom ligands here (50% regular passivation)
    for i in range(0, nx, 2):
        for j in range(0, nz):
            coords = np.copy(centred_bot)
            # shift molecule
            coords[:, 0] = coords[:, 0] + 4.12 * i
            coords[:, 1] = coords[:, 1] + -15.457550
            coords[:, 2] = coords[:, 2] + 6.75 * j + 3.375
#             print site_x, site_z        

            for k in range(n_chain):
                if (k == n_chain-1):
                    print "4       {:4.5f} {:4.5f} {:4.5f}".format( coords[k,0], coords[k,1], coords[k,2] )
                elif ((k == 7) or (k == 10)):
                    print "2       {:4.5f} {:4.5f} {:4.5f}".format( coords[k,0], coords[k,1], coords[k,2] )
                elif ((k == 8) or (k == 9)):
                    print "3       {:4.5f} {:4.5f} {:4.5f}".format( coords[k,0], coords[k,1], coords[k,2] )
                else:
                    print "1       {:4.5f} {:4.5f} {:4.5f}".format( coords[k,0], coords[k,1], coords[k,2] )
    
            nmol = nmol + 1

    # hex molecules here

    coords = centred_hex + np.array([0.0, 22.0, 0.0])
    for i in range(n_hexatoms):    
        if ((i%6 == 0) or ((i+1)%6 == 0)):
            print "7	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[i,0], coords[i,1], coords[i,2] )
        else:
            print "8	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[i,0], coords[i,1], coords[i,2] )

if args.fmt == "lmp":
    # beginning of LAMMPS input file
    print "LAMMPS coords for {} oleic acid-like molecules\n".format( n_molecules )
    print "{}	atoms".format( n_top * n_chain + n_bottom * n_chain + n_molecules * 16 + n_hexatoms )
    print "{}	bonds".format( (n_top + n_bottom) * n_bonds + n_top + n_bottom + n_hexmolecules * n_hexbonds )
    print "{}	angles".format( (n_top + n_bottom) * n_angles + n_hexmolecules * n_hexangles )
    print "{}	dihedrals".format( (n_top + n_bottom) * n_dihedrals + n_hexmolecules * n_hexdihedrals )
    print "{}	impropers\n".format( n_top + n_bottom )
    
    print "8	atom types"
    print "2	bond types"
    print "2	angle types"
    print "2	dihedral types"
    print "1	improper types\n"
    
    print "{} {}	xlo xhi".format(xlo, xhi)
    print "{} {}	ylo yhi".format(ylo, yhi)
    print "{} {}	zlo zhi\n".format(zlo, zhi)
    
    print "Masses\n"
    print "1 14.027"
    print "2 14.027"
    print "3 13.019"
    print "4 15.034"
    print "5 112.411"
    print "6 32.065"
    print "7 15.034"
    print "8 14.027\n"
    print "Atoms\n"

    num_mol = 1
    num_atom = 1
    N_lig = int(args.cov * nx * nz)
    occ_indices = []
    while (N_lig > 0):
        site_x = randint(0, nx - 1)
        site_z = randint(0, nz - 1)
        site_num = [site_x, site_z]
    
        if site_num not in occ_indices:
            # place ligand
            coords = np.copy(centred_data)
            # shift molecule
            coords[:, 0] = coords[:, 0] + 4.12 * site_x
            coords[:, 1] = coords[:, 1]
            coords[:, 2] = coords[:, 2] + 6.75 * site_z
    
            for k in range(n_chain):
                if (k == n_chain-1):
                    print "{}	{}	4       {:4.5f} {:4.5f} {:4.5f}".format( num_atom, num_mol, coords[k,0], coords[k,1], coords[k,2] )
                elif ((k == 7) or (k == 10)):
                    print "{}	{}	2       {:4.5f} {:4.5f} {:4.5f}".format( num_atom, num_mol, coords[k,0], coords[k,1], coords[k,2] )
                elif ((k == 8) or (k == 9)):
                    print "{}	{}	3       {:4.5f} {:4.5f} {:4.5f}".format( num_atom, num_mol, coords[k,0], coords[k,1], coords[k,2] )
                else:
                    print "{}	{}	1       {:4.5f} {:4.5f} {:4.5f}".format( num_atom, num_mol, coords[k,0], coords[k,1], coords[k,2] )
                num_atom = num_atom + 1

                if ((k+1)%n_chain == 0):
                    num_mol = num_mol + 1

            occ_indices.append(site_num)
#             mask[site_x, site_z] = 1
            N_lig = N_lig - 1
#     occ_indices = np.array(occ_indices)

    # Cd-S surface coordinates
    Cd_top = np.zeros((nx, nz))
    Cd_bot = []
    for i in range(nx):
        for j in range(nz):
            coords = centred_data[0, :] - np.array([0.0, 1.56, 0.0])
            coords[0] = coords[0] + 4.12 * i
            coords[2] = coords[2] + 6.75 * j
            print "{}	{}	5	{:4.5f}	{:4.5f}	{:4.5f}".format( num_atom, num_mol, coords[0], coords[1], coords[2] )
            Cd_top[i, j] = num_atom
            num_atom = num_atom + 1
            print "{}	{}	6	{:4.5f}	{:4.5f}	{:4.5f}".format( num_atom, num_mol, coords[0], coords[1], coords[2] + 2.5299 )
            num_atom = num_atom + 1

    for i in range(nx):
        for j in range(nz):
            coords = centred_data[0, :] - np.array([0.0, 1.56, 0.0])
            coords[0] = coords[0] + 4.12 * i
            coords[2] = coords[2] + 6.75 * j
            print "{}	{}	5	{:4.5f}	{:4.5f}	{:4.5f}".format( num_atom, num_mol, coords[0] - 2.06073, coords[1] - 1.18977, coords[2] + 3.375 )
            num_atom = num_atom + 1
            print "{}	{}	6	{:4.5f}	{:4.5f}	{:4.5f}".format( num_atom, num_mol, coords[0] - 2.06073, coords[1] - 1.18977, coords[2] - 3.375 + 2.5299 )
            num_atom = num_atom + 1
            print "{}	{}	5	{:4.5f}	{:4.5f}	{:4.5f}".format( num_atom, num_mol, coords[0] - 2.06073, coords[1] - 3.5693, coords[2] )
            num_atom = num_atom + 1
            print "{}	{}	6	{:4.5f}	{:4.5f}	{:4.5f}".format( num_atom, num_mol, coords[0] - 2.06073, coords[1] - 3.5693, coords[2] + 2.5299 )
            num_atom = num_atom + 1
            print "{}	{}	5	{:4.5f}	{:4.5f}	{:4.5f}".format( num_atom, num_mol, coords[0], coords[1] - 3.5693 - 1.18977, coords[2] + 3.375 )
            num_atom = num_atom + 1
            print "{}	{}	6	{:4.5f}	{:4.5f}	{:4.5f}".format( num_atom, num_mol, coords[0], coords[1] - 3.5693 - 1.18977, coords[2] - 3.375 + 2.5299 )
            num_atom = num_atom + 1

            print "{}	{}	5	{:4.5f}	{:4.5f}	{:4.5f}".format( num_atom, num_mol, coords[0], coords[1] - 7.13848, coords[2] )
            num_atom = num_atom + 1
            print "{}	{}	6	{:4.5f}	{:4.5f}	{:4.5f}".format( num_atom, num_mol, coords[0], coords[1] - 7.13848, coords[2] + 2.5299 )
            num_atom = num_atom + 1
            print "{}	{}	5	{:4.5f}	{:4.5f}	{:4.5f}".format( num_atom, num_mol, coords[0] - 2.06073, coords[1] - 7.13848 - 1.18977, coords[2] + 3.375 )
            num_atom = num_atom + 1
            print "{}	{}	6	{:4.5f}	{:4.5f}	{:4.5f}".format( num_atom, num_mol, coords[0] - 2.06073, coords[1] - 7.13848 - 1.18977, coords[2] - 3.375 + 2.5299 )
            num_atom = num_atom + 1
            print "{}	{}	5	{:4.5f}	{:4.5f}	{:4.5f}".format( num_atom, num_mol, coords[0] - 2.06073, coords[1] - 7.13848 - 3.5693, coords[2] )
            num_atom = num_atom + 1
            print "{}	{}	6	{:4.5f}	{:4.5f}	{:4.5f}".format( num_atom, num_mol, coords[0] - 2.06073, coords[1] - 7.13848 - 3.5693, coords[2] + 2.5299 )
            num_atom = num_atom + 1
            print "{}	{}	5	{:4.5f}	{:4.5f}	{:4.5f}".format( num_atom, num_mol, coords[0], coords[1] - 7.13848 - 3.5693 - 1.18977, coords[2] + 3.375 )
            Cd_bot.append(num_atom)
            num_atom = num_atom + 1
            print "{}	{}	6	{:4.5f}	{:4.5f}	{:4.5f}".format( num_atom, num_mol, coords[0], coords[1] - 7.13848 - 3.5693 - 1.18977, coords[2] - 3.375 + 2.5299 )
            num_atom = num_atom + 1

    Cd_bot = np.array(Cd_bot)
#     print "bottom indices"
#     print Cd_bot

    # bottom ligands at 50% passivation
    mask = np.zeros((nx, nz))
    for i in range(nx):
        for j in range(0, nz, 2):
            coords = np.copy(centred_bot)
            # shift molecule
            mask[i, j] = 1
            coords[:, 0] = coords[:, 0] + 4.12 * i
            coords[:, 1] = coords[:, 1] + -15.457550
            coords[:, 2] = coords[:, 2] + 6.75 * j + 3.375
        
            for k in range(n_chain):
                if (k == n_chain-1):
                    print "{}	{}	4	{:4.5f}	{:4.5f}	{:4.5f}".format( num_atom, num_mol, coords[k,0], coords[k,1], coords[k,2]  )
                elif ((k == 7) or (k == 10)):
                    print "{}	{}	2	{:4.5f}	{:4.5f}	{:4.5f}".format( num_atom, num_mol, coords[k,0], coords[k,1], coords[k,2]  )
                elif ((k == 8) or (k == 9)):
                    print "{}	{}	3	{:4.5f}	{:4.5f}	{:4.5f}".format( num_atom, num_mol, coords[k,0], coords[k,1], coords[k,2]  )
                else:
                    print "{}	{}	1	{:4.5f}	{:4.5f}	{:4.5f}".format( num_atom, num_mol, coords[k,0], coords[k,1], coords[k,2]  )
                num_atom = num_atom + 1 

            num_mol = num_mol + 1

#     print "bottom ligand mask"
#     print mask

    coords = centred_hex + np.array([0.0, 22.0, 0.0]) 
    for i in range(n_hexatoms):
        if ((i%6 == 0) or ((i+1)%6 == 0)):
            print "{}	{}	7	{:4.5f}	{:4.5f}	{:4.5f}".format( num_atom, num_mol, coords[i,0], coords[i,1], coords[i,2] )
        else:
            print "{}	{}	8	{:4.5f}	{:4.5f}	{:4.5f}".format( num_atom, num_mol, coords[i,0], coords[i,1], coords[i,2] )

        num_atom = num_atom + 1
        # increase molecule number every 6 atoms
        if ((i+1)%6 == 0):
            num_mol = num_mol + 1

    # start defining how atoms are connected
    print "\nBonds\n"

    # define C-C and C=C bonds here for top ligands
    num_bond = 1
    for i in range(n_top):
        for j in range(n_bonds):
            first = i * n_chain + j + 1
            if j == 8:
                print "{}	2	{}	{}".format( num_bond, first, first + 1 )
            else:
                print "{}	1	{}	{}".format( num_bond, first, first + 1 )
            num_bond = num_bond + 1

#     print "finished upper bonds" 
    # define Cd-CH2 bonds here for top ligands
    for i in range(n_top):
        carbon = i * n_chain + 1
        site_num = occ_indices[i]
        cadmium = int( Cd_top[site_num[0], site_num[1]] )
        print "{}	1	{}	{}".format( num_bond, carbon, cadmium )
        num_bond = num_bond + 1

#     print "finished upper cd-c bonds" 
    # define C-C and C=C bonds here for bottom ligands
    for i in range(n_bottom):
        for j in range(n_bonds):
            first = n_top * n_chain + n_molecules * 16 + i * n_chain + j + 1
            if j == 8:
                print "{}	2	{}	{}".format( num_bond, first, first + 1 )
            else:
                print "{}	1	{}	{}".format( num_bond, first, first + 1 )
            num_bond = num_bond + 1

#     print "finished lower bonds" 
    # define Cd-CH2 bonds here for top ligands
    for i in range(n_bottom):
        carbon = n_top * n_chain + n_molecules * 16 + i * n_chain + 1
        cadmium = Cd_bot[2 * i] 
        print "{}	1	{}	{}".format( num_bond, carbon, cadmium )
        num_bond = num_bond + 1

    # define hex bonds here
    for i in range(n_hexmolecules):
        for j in range(n_hexbonds):
            first = (n_top + n_bottom) * n_chain + n_molecules * 16 + i * n_hex + j + 1
            print "{}	1	{}	{}".format( num_bond, first, first + 1 )
            num_bond = num_bond + 1

    print "\nAngles\n"

    num_angle = 1
    # top oleic acid angles here
    for i in range(n_top):
        for j in range(n_angles):
            first = i * n_chain + j + 1
            if (j == 7) or (j == 8):
                print "{}	2	{}	{}	{}".format( num_angle, first, first + 1, first + 2 )
            else:
                print "{}	1	{}	{}	{}".format( num_angle, first, first + 1, first + 2 )
            num_angle = num_angle + 1

    # bottom oleic acid angles here
    for i in range(n_bottom):
        for j in range(n_angles):
            first = n_top * n_chain + n_molecules * 16 + i * n_chain + j + 1
            if (j == 7) or (j == 8):
                print "{}	2	{}	{}	{}".format( num_angle, first, first + 1, first + 2 )
            else:
                print "{}	1	{}	{}	{}".format( num_angle, first, first + 1, first + 2 )
            num_angle = num_angle + 1


    # hex angles defined here
    for i in range(n_hexmolecules):
        for j in range(n_hexangles):
            first = (n_top + n_bottom) * n_chain + n_molecules * 16 + i * n_hex + j + 1
            print "{}	1	{}	{}	{}".format( num_angle, first, first + 1, first + 2)
            num_angle = num_angle + 1
    
    print "\nDihedrals\n"

    num_dihedral = 1
    # top oleic acid dihedrals here
    for i in range(n_top):
        for j in range(n_dihedrals):
            if j > 6:
                first = i * n_chain + j + 2
            else:
                first = i * n_chain + j + 1
    
            if (j == 6) or (j == 7):
                print "{}	2	{}	{}	{}	{}".format( num_dihedral, first, first + 1, first + 2, first + 3 )
            else:
                print "{}	1	{}	{}	{}	{}".format( num_dihedral, first, first + 1, first + 2, first + 3 )
            num_dihedral = num_dihedral + 1

    # bottom oleic acid dihedrals here
    for i in range(n_bottom):
        for j in range(n_dihedrals):
            if j > 6:
                first = n_top * n_chain + n_molecules * 16 + i * n_chain + j + 2
            else:
                first = n_top * n_chain + n_molecules * 16 + i * n_chain + j + 1
    
            if (j == 6) or (j == 7):
                print "{}	2	{}	{}	{}	{}".format( num_dihedral, first, first + 1, first + 2, first + 3 )
            else:
                print "{}	1	{}	{}	{}	{}".format( num_dihedral, first, first + 1, first + 2, first + 3 )
            num_dihedral = num_dihedral + 1


    # hex dihedrals here
    for i in range(n_hexmolecules):
        for j in range(n_hexdihedrals):
            first = (n_top + n_bottom) * n_chain + n_molecules * 16 + i * n_hex + j + 1
            print "{}	1	{}	{}	{}	{}".format( num_dihedral, first, first + 1, first + 2, first + 3 )
            num_dihedral = num_dihedral + 1

    print "\nImpropers\n"
    
    for i in range(n_top):
        first = i * n_chain + 8
        print "{}	1	{}	{}	{}	{}".format( i + 1, first, first + 1, first + 2, first + 3 )

    for i in range(n_bottom):
        first = n_top * n_chain + n_molecules * 16 + i * n_chain + 8
        print "{}	1	{}	{}	{}	{}".format( n_top + i + 1, first, first + 1, first + 2, first + 3 )


