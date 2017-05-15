from __future__ import print_function
import numpy as np

log = np.genfromtxt('log.stripped').astype(int)

filename = ["pot-349.5.txt",
	    "pot-350.0.txt",
	    "pot-350.5.txt",
	    "pot-351.0.txt",
	    "pot-351.5.txt",
	    "pot-352.0.txt",
            "pot-352.5.txt",
            "pot-353.0.txt"]

newfiles = ["pot-new.349.5.txt",
	    "pot-new.350.0.txt",
	    "pot-new.350.5.txt",
	    "pot-new.351.0.txt",
	    "pot-new.351.5.txt",
	    "pot-new.352.0.txt",
            "pot-new.352.5.txt",
            "pot-new.353.0.txt"]

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

