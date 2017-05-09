from __future__ import print_function
import numpy as np

log = np.genfromtxt('log.stripped').astype(int)

filename = ["dump.0.349.0",
	    "dump.0.349.5",
	    "dump.0.350.0",
	    "dump.0.350.5",
	    "dump.0.351.0",
	    "dump.0.351.5"]

newfiles = ["dump.new.349.0",
	    "dump.new.349.5",
	    "dump.new.350.0",
	    "dump.new.350.5",
	    "dump.new.351.0",
	    "dump.new.351.5"]

openfiles = []

for f in filename:
    openfiles.append(open(f, 'r'))

first_time = True
for row in log:
    time = row[0]
    if (time%10000 != 0):
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

