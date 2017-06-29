from __future__ import print_function
import numpy as np

log = np.genfromtxt('log.stripped').astype(int)
t_max = [None, None, None, None, None, None, None, None]
tdist = []

for row in log:
    time = row[0]
#     print(time)
    minint = row[1]
    maxint = row[8]
    if t_max[maxint] is None:
        t_max[maxint] = time
    if not t_max[minint] is None:
        tdist.append(time - t_max[minint])
        t_max[minint] = None

print(t_max)
print(time)

