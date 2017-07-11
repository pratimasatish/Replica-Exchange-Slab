import numpy as np
import scipy.special as sp
import scipy.optimize as opt
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("file_in", type=str, help="Name of file to load data from")
parser.add_argument('-v', action='store_true', help="Run verbosely")
args = parser.parse_args()


def bessel_0(x):
    return sp.kn(0, x)

def chi_func(F, L, x, y, dy):
    assert( len(x) == len(y) )
    chisq = 0
    for i in xrange(len(x)):
        err   = y[i] - F(x[i] / float(L) )
        err_s = err / dy[i]
        chisq += err_s * err_s
    chisq /= len(x)
    return chisq

sheet_data = np.loadtxt(args.file_in)
x_max = max(sheet_data[:,0])
z_max = max(sheet_data[:,1])
corr_xz = np.zeros((x_max+1, z_max+1))
dcorr_xz = np.zeros((x_max+1, z_max+1))
for row in sheet_data:
    corr_xz[row[0], row[1]] = row[2]
    dcorr_xz[row[0], row[1]] = row[3]

# choose either x- or z-direction for fitting

#axis = range(1,int(x_max+1))  # x-direction
#dat  = corr_xz[1:, 0]
#d_dat = dcorr_xz[1:, 0]

axis = range(1,int(z_max+1))   # z-direction
dat  = corr_xz[0, 1:]
d_dat = dcorr_xz[0, 1:]

def chi_bessel(L):
    return chi_func(bessel_0, L, axis, dat, d_dat)
# optL = opt.minimize(chi_bessel, .1, method='TNC', bounds=((0,None),) ).x[0]
optL = opt.minimize(chi_bessel, .1).x[0]

def chi_shift(L):
    return abs(chi_bessel(L) - chi_bessel(optL) - 1)

# L_test = np.linspace(optL, 2*optL, 19999)
# chi_test = [chi_shift(Lt) for Lt in L_test]
# plt.plot(L_test, chi_test)
# plt.show()

varL = opt.minimize(chi_shift, optL*1.1, method='TNC', bounds=((optL,None),) ).x[0]
errL = varL - optL
# print varL, errL

if args.v:
    print "Minimum value of L +/- dL for K_0(x/L) fit to data:",
print "{} {}".format( optL, errL )
