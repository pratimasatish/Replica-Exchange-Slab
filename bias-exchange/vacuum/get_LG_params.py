import numpy as np

# in K
temp = np.array([370, 380])

# in rad
ord_min = np.array([-0.727, -0.714])
disord_min = np.array([-0.326, -0.287])

# in kcal/mol-K
del_S = np.array([1.313, 1.418])

d0 = -3.0 * del_S / (disord_min**2 - ord_min**2)
avg_d0 = np.mean(d0)

print "d0 at 370 K: {} kcal/mol-K-rad^2".format(d0[0])
print "d0 at 380 K: {} kcal/mol-K-rad^2".format(d0[1])
print "average d0: {} kcal/mol-K-rad^2".format(avg_d0)

# in kcal/mol-K-nm
sigma = 0.0129

# in nm
lambda_val = np.array([9.72, 10.8])
lambda_sq = lambda_val ** 2
eta = np.array([-0.73, -0.24])
T = np.array([365, 400])

 



