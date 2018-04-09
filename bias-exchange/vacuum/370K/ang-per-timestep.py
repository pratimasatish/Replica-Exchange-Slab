import numpy as np
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("-bias", type=str, help="bias value to analyse")
args = parser.parse_args()

data = np.genfromtxt('theta-ag' + args.bias + '.txt')
f = open("theta-ag-time" + args.bias + ".txt","w")

data = np.mean(data.reshape((-1, 240)), axis=1)

for i in range(len(data)):
    f.write("%4.5f\n" %(data[i]))

f.close()

