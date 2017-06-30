from __future__ import print_function
import numpy as np

log = np.genfromtxt('log.stripped').astype(int)

filename = ["en-340.txt",
	    "en-343.txt",
            "en-346.txt",
	    "en-349.txt",
	    "en-352.txt",
	    "en-355.txt",
	    "en-358.txt",
            "en-361.txt",
            "en-364.txt",
            "en-367.txt",
            "en-370.txt"]

newfiles = ["en-new.340.txt",
	    "en-new.343.txt",
	    "en-new.346.txt",
	    "en-new.349.txt",
	    "en-new.352.txt",
	    "en-new.355.txt",
            "en-new.358.txt",
            "en-new.361.txt",
            "en-new.364.txt",
            "en-new.367.txt",
            "en-new.370.txt"]

openfiles = []

for f in filename:
    openfiles.append(open(f, 'r'))

first_time = True
for row in log:
    time = row[0] - 1
    if (time%10000 != 0):
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

