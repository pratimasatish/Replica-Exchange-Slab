from __future__ import print_function
import numpy as np


log = np.genfromtxt('log.stripped').astype(int)

filename = ["dump.temper.340.0",
	    "dump.temper.343.0",
	    "dump.temper.346.0",
	    "dump.temper.349.0",
	    "dump.temper.352.0",
	    "dump.temper.355.0",
	    "dump.temper.358.0",
	    "dump.temper.361.0",
	    "dump.temper.364.0",
	    "dump.temper.367.0",
	    "dump.temper.370.0"]

newfiles = ["dump.new.340.xyz",
	    "dump.new.343.xyz",
	    "dump.new.346.xyz",
	    "dump.new.349.xyz",
	    "dump.new.352.xyz",
	    "dump.new.355.xyz",
	    "dump.new.358.xyz",
	    "dump.new.361.xyz",
	    "dump.new.364.xyz",
	    "dump.new.367.xyz",
	    "dump.new.370.xyz"]

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

