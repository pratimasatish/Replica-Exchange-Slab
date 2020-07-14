import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("-temp", default='365', type=str, help="temp value to analyse")
args = parser.parse_args()

plt.rc('text', usetex=True)
plt.rc('font', family='serif', weight='bold')
plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
plt.rc('font', size=28)
plt.rc('xtick', labelsize=28)
plt.rc('ytick', labelsize=28)

data_z = np.genfromtxt('theta' + args.temp + '.txt')
data_x = np.genfromtxt('theta-x' + args.temp + '.txt')

bin_z = np.linspace(data_z.min(), data_z.max(), 100)
bin_x = np.linspace(data_x.min(), data_x.max(), 100)

hist = np.histogram2d(data_z, data_x, bins=(bin_z, bin_x), density=True)[0]
hist = hist.T
zz, xx = np.meshgrid(bin_z, bin_x)
print(hist.shape)

plt.pcolormesh(zz, xx, hist, cmap="Reds")
plt.axis( [zz.min(), zz.max(), xx.min(), xx.max()] )
plt.colorbar()
plt.title(r'$\textrm{Joint probability of } \theta_z \textrm{ and } \theta_x$')
plt.show()

inds = np.array( [0, 12, 26, 33, 48, 55, 64, 79, 87, 95] )
# chosen_z = np.array( [bin_z[0], bin_z[12], bin_z[26], bin_z[33], bin_z[48], bin_z[55], bin_z[64], bin_z[79], bin_z[87], bin_z[95]] )
for i in range(len(inds)):
    print('at theta_z = {} , the mean is {} and the std dev is {}').format( bin_z[inds[i]], np.mean(hist[inds[i], :]), np.std(hist[inds[i], :]) )


data_tz = np.reshape(data_z, (-1, 240))
ave_z = np.mean(data_tz, axis=1)
data_tx = np.reshape(data_x, (-1, 240))
ave_x = np.mean(data_tx, axis=1)

print ave_x.shape
print('<theta_z> mean = {} and standard dev = {}').format(np.mean(ave_x), np.std(ave_x))

avebin_z = np.linspace(ave_z.min(), ave_z.max(), 100)
avebin_x = np.linspace(ave_x.min(), ave_x.max(), 100)

avehist = np.histogram2d(ave_z, ave_x, bins=(avebin_z, avebin_x), density=True)[0]
avehist = avehist.T
ave_zz, ave_xx = np.meshgrid(avebin_z, avebin_x)


plt.pcolormesh(ave_zz, ave_xx, avehist, cmap="Reds")
# plt.imshow(avehist, cmap="Reds", origin="lower", interpolation="none")
plt.axis( [ave_zz.min(), ave_zz.max(), ave_xx.min(), ave_xx.max()] )
plt.colorbar()
plt.title(r'$\textrm{Joint probability of } \langle\theta_z\rangle \textrm{ and } \langle\theta_x\rangle$')
plt.show()

hist_x = np.histogram(data_x, bins=bin_x, density=True)[0]
bin_centres = 0.5 * (bin_x[1:] + bin_x[:-1])
plt.plot(bin_centres, hist_x, 'ro')
plt.plot(bin_centres, hist_x, 'r', linewidth=4, alpha=0.7)
plt.title(r'$P(\theta_x)$')
plt.show()

avehist_x = np.histogram(ave_x, bins=avebin_x, density=True)[0]
ave_centres = 0.5 * (avebin_x[1:] + avebin_x[:-1])
plt.plot(ave_centres, avehist_x, 'ro')
plt.plot(ave_centres, avehist_x, 'r', linewidth=4, alpha=0.7)
plt.title(r'$P(\langle\theta_x\rangle)$')
plt.show()

# data_ord_z = np.genfromtxt('theta-0.750.txt')
# data_ord_x = np.genfromtxt('theta-x-0.750.txt')
# data_disord_z = np.genfromtxt('theta-0.290.txt')
# data_disord_x = np.genfromtxt('theta-x-0.290.txt')
# 
# bin_ord_x = np.linspace(data_ord_x.min(), data_ord_x.max(), 100)
# bin_disord_x = np.linspace(data_disord_x.min(), data_disord_x.max(), 100)
# 
# hist_ord_x = np.histogram(data_ord_x, bins=bin_ord_x, density=True)[0]
# bin_ord_centres = 0.5 * (bin_ord_x[1:] + bin_ord_x[:-1])
# hist_disord_x = np.histogram(data_disord_x, bins=bin_disord_x, density=True)[0]
# bin_disord_centres = 0.5 * (bin_disord_x[1:] + bin_disord_x[:-1])
# plt.plot(bin_ord_centres, hist_ord_x, 'ro')
# plt.plot(bin_ord_centres, hist_ord_x, 'r', linewidth=4, alpha=0.7, label='ordered')
# plt.plot(bin_disord_centres, hist_disord_x, 'bv')
# plt.plot(bin_disord_centres, hist_disord_x, 'b', linewidth=4, alpha=0.7, label='disordered')
# plt.legend(loc='best')
# plt.title(r'$P(\theta_x)$')
# plt.show()
# 
# data_ord_tx = np.reshape(data_ord_x, (-1, 240))
# ave_ord_x = np.mean(data_ord_tx, axis=1)
# data_disord_tx = np.reshape(data_disord_x, (-1, 240))
# ave_disord_x = np.mean(data_disord_tx, axis=1)
# 
# avebin_ord_x = np.linspace(ave_ord_x.min(), ave_ord_x.max(), 100)
# avebin_disord_x = np.linspace(ave_disord_x.min(), ave_disord_x.max(), 100)
# 
# avehist_ord_x = np.histogram(ave_ord_x, bins=avebin_ord_x, density=True)[0]
# avehist_disord_x = np.histogram(ave_disord_x, bins=avebin_disord_x, density=True)[0]
# ave_ord_centres = 0.5 * (avebin_ord_x[1:] + avebin_ord_x[:-1])
# ave_disord_centres = 0.5 * (avebin_disord_x[1:] + avebin_disord_x[:-1])
# plt.plot(ave_ord_centres, avehist_ord_x, 'ro', label=r'$P(\langle\theta_x\rangle)_{ord}$')
# plt.plot(ave_ord_centres, avehist_ord_x, 'r', linewidth=4, alpha=0.7)
# plt.plot(ave_disord_centres, avehist_disord_x, 'bv', label=r'$P(\langle\theta_x\rangle)_{disord}$')
# plt.plot(ave_disord_centres, avehist_disord_x, 'b', linewidth=4, alpha=0.7)
# plt.legend(loc='best')
# plt.title(r'$P(\langle\theta_x\rangle)$')
# plt.show()
# 
