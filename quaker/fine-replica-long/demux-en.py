from __future__ import print_function
import numpy as np

log = np.genfromtxt('log.stripped').astype(int)

filename = ["en-349.0.txt",
	    "en-349.5.txt",
	    "en-350.0.txt",
	    "en-350.5.txt",
	    "en-351.0.txt",
	    "en-351.5.txt"]

newfiles = ["en-new.349.0.txt",
	    "en-new.349.5.txt",
	    "en-new.350.0.txt",
	    "en-new.350.5.txt",
	    "en-new.351.0.txt",
	    "en-new.351.5.txt"]

openfiles = []

for f in filename:
    openfiles.append(open(f, 'r'))

first_time = True
for row in log:
    time = row[0]
    if (time%1000 != 0):
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

