from __future__ import print_function
import numpy as np


log = np.genfromtxt('log.stripped').astype(int)

filename = ["dump.-0.750.xyz",
	    "dump.-0.700.xyz",
	    "dump.-0.650.xyz",
	    "dump.-0.625.xyz",
	    "dump.-0.600.xyz",
	    "dump.-0.580.xyz",
	    "dump.-0.560.xyz",
	    "dump.-0.540.xyz",
	    "dump.-0.520.xyz",
	    "dump.-0.500.xyz",
	    "dump.-0.425.xyz",
	    "dump.-0.350.xyz",
	    "dump.-0.275.xyz",
	    "dump.-0.225.xyz",
	    "dump.-0.175.xyz",
            "dump.-0.125.xyz"]

newfiles = ["dump.new.-0.750.xyz",
	    "dump.new.-0.700.xyz",
	    "dump.new.-0.650.xyz",
	    "dump.new.-0.625.xyz",
	    "dump.new.-0.600.xyz",
	    "dump.new.-0.580.xyz",
	    "dump.new.-0.560.xyz",
	    "dump.new.-0.540.xyz",
	    "dump.new.-0.520.xyz",
	    "dump.new.-0.500.xyz",
	    "dump.new.-0.425.xyz",
	    "dump.new.-0.350.xyz",
	    "dump.new.-0.275.xyz",
	    "dump.new.-0.225.xyz",
	    "dump.new.-0.175.xyz",
            "dump.new.-0.125.xyz"]

openfiles = []

for f in filename:
    openfiles.append(open(f, 'r'))

for f in openfiles:
    [f.readline() for frame_count in range(10322)]

first_time = True
# for row in log:
for j in range(0,len(log)-1):
    row = log[j]
    time = row[0]
#     print(j, row, time)
    if (time%5000 != 0):
#     if (time%1000 != 0):
        continue
    print(time)
    for i,f in enumerate(openfiles):
        print(i, row[i+1])
        lines = [f.readline() for frame_count in range(10322)]
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

