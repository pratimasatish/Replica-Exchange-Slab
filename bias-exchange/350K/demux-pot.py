from __future__ import print_function
import numpy as np

log = np.genfromtxt('log.lammps').astype(int)

filename = ["pot-0.750",
	    "pot-0.650",
	    "pot-0.600",
	    "pot-0.575",
            "pot-0.550",
            "pot-0.525",
            "pot-0.500",
            "pot-0.450",
            "pot-0.375",
            "pot-0.275",
            "pot-0.200",
            "pot-0.100"]

newfiles = ["pot-new.-0.750",
	    "pot-new.-0.650",
	    "pot-new.-0.600",
	    "pot-new.-0.575",
            "pot-new.-0.550",
            "pot-new.-0.525",
            "pot-new.-0.500",
            "pot-new.-0.450",
            "pot-new.-0.375",
            "pot-new.-0.275",
            "pot-new.-0.200",
            "pot-new.-0.100"]

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

