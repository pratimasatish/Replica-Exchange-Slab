import numpy as np
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="")
parser.add_argument("-temp", type=str, help="temp value to analyse")
parser.add_argument("-step", type=int, default=500, help="step value for running mean")
args = parser.parse_args()

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / N

data = np.genfromtxt('theta' + args.temp + '.txt', delimiter=' ')
data = np.mean( data.reshape((-1, 240)), axis=1 )

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.figure()
plt.plot(running_mean(data, args.step))
plt.ylim(-0.85, -0.25)
plt.show()

