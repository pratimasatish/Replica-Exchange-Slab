import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-nx", type=int, default=2, help="number of ligands in x-direction")
parser.add_argument("-nz", type=int, default=2, help="number of ligands in z-direction")
parser.add_argument("-fmt", type=str, choices=["xyz", "lmp"], default='xyz', help="output format")
args = parser.parse_args()

data = np.genfromtxt('vacuum-frame.txt')

# number of replicas to be made
nreps = args.nx * args.nz
n_atoms = len(data)

Lx = 82.4293
Lz = 81.004
xlo =  0.0
xhi =  Lx * args.nx
ylo = -60.0
yhi = 120.0
zlo =  0.0
zhi =  Lz * args.nz
n_top = 240
n_bottom = 120
n_chain = 18
n_bonds = n_chain - 1
n_angles = n_chain - 2
n_dihedrals = n_chain - 3

if args.fmt == 'xyz':
    # start replicating the box
    print "{}".format(n_atoms * nreps)
    print "Comment Line"

    for i in range(args.nx):
        for j in range(args.nz):
            for k in range(n_atoms):
                print "{} {:4.5f} {:4.5f} {:4.5f}".format( data[k, 0], data[k, 1] + i * Lx, data[k, 2], data[k, 3] + j * Lz )

if args.fmt == 'lmp':
    print "LAMMPS coords for {} extended octadecyl slabs\n".format( nreps )
    print "{}	atoms".format( n_atoms * nreps )
    print "{}	bonds".format( ( (n_top + n_bottom) * n_bonds + n_top + n_bottom ) * nreps )
    print "{}	angles".format( ((n_top + n_bottom) * n_angles) * nreps )
    print "{}	dihedrals".format( ((n_top + n_bottom) * n_dihedrals) * nreps )
    print "0	impropers\n"
    
    print "6	atom types"
    print "1	bond types"
    print "1	angle types"
    print "1	dihedral types"
    print "0	improper types\n"
    
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


    exit(1)

