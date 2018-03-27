from __future__ import print_function
import numpy as np

log = np.genfromtxt('log.lammps').astype(int)

filename = ["pot.-0.750",
	    "pot.-0.700",
	    "pot.-0.650",
	    "pot.-0.625",
	    "pot.-0.600",
	    "pot.-0.580",
            "pot.-0.560",
            "pot.-0.540",
            "pot.-0.520",
            "pot.-0.500",
            "pot.-0.425",
            "pot.-0.350",
            "pot.-0.275",
            "pot.-0.225",
            "pot.-0.175",
            "pot.-0.125"]

newfiles = ["pot-new.-0.750",
	    "pot-new.-0.700",
	    "pot-new.-0.650",
	    "pot-new.-0.625",
	    "pot-new.-0.600",
	    "pot-new.-0.580",
            "pot-new.-0.560",
            "pot-new.-0.540",
            "pot-new.-0.520",
            "pot-new.-0.500",
            "pot-new.-0.425",
            "pot-new.-0.350",
            "pot-new.-0.275",
            "pot-new.-0.225",
            "pot-new.-0.175",
            "pot-new.-0.125"]

openfiles = []

for f in filename:
    openfiles.append(open(f, 'r'))

first_time = True
for row in log:
    time = row[0]
    if (time%5000 != 0):
        continue
    print(time)
    for i,f in enumerate(openfiles):
        print(i, row[i+1])
        lines = f.readline()
        if first_time:
            g = open(newfiles[row[i+1]], 'w+')
        else:
            g = open(newfiles[row[i+1]], 'a')
	g.write(lines)
        g.close()
    first_time=False

for f in openfiles:
    f.close()

