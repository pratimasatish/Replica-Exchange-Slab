from __future__ import print_function
import numpy as np


log = np.genfromtxt('log.stripped').astype(int)

filename = ["dump.rplc_plm.355.0.xyz",
	    "dump.rplc_plm.358.0.xyz"]

newfiles = ["dump.new.355.xyz",
	    "dump.new.358.xyz"]

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
        lines = [f.readline() for frame_count in range(22118)]
	lines = ''.join(lines)
        if first_time:
            g = open(newfiles[row[i+1]], 'w+')
        else:
            g = open(newfiles[row[i+1]], 'a')
	g.write(lines)
        g.close()
    first_time=False

for f in openfiles:
    f.close()

