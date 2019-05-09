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

plt.pcolormesh(zz, xx, hist, cmap="Reds")
plt.axis( [zz.min(), zz.max(), xx.min(), xx.max()] )
plt.colorbar()
plt.title(r'$\textrm{Joint probability of } \theta_z \textrm{ and } \theta_x$')
plt.show()

data_tz = np.reshape(data_z, (-1, 240))
ave_z = np.mean(data_tz, axis=1)
data_tx = np.reshape(data_x, (-1, 240))
ave_x = np.mean(data_tx, axis=1)

avebin_z = np.linspace(ave_z.min(), ave_z.max(), 100)
avebin_x = np.linspace(ave_x.min(), ave_x.max(), 100)

avehist = np.histogram2d(ave_z, ave_x, bins=(avebin_z, avebin_x), density=True)[0]
avehist = avehist.T
ave_zz, ave_xx = np.meshgrid(avebin_z, avebin_x)

plt.pcolormesh(ave_zz, ave_xx, avehist, cmap="Reds")
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


