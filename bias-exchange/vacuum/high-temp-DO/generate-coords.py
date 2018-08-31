import numpy as np

T = 500
n_particle = 100
L_box = 10

coords = L_box * np.random.random_sample((T * n_particle, 3))
# print coords.shape
# print len(coords)

for i in range(len(coords)):
    print "1 {} {} {}".format(coords[i][0], coords[i][1], coords[i][2])
