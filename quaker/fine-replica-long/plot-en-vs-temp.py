import numpy as np
import matplotlib.pyplot as plt

def read_file(filename):
    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()
    return lines

temperature_k = [349.0, 349.5, 350.0, 350.5, 351.0, 351.5]
K = len(temperature_k)
T = 501

U_kn = np.zeros((K, T))		# U_kn[k,t] is the potential energy (in kcal/mol) for snapshot t of temperature index k
for k in range(K):
    lines = read_file("en-new.{:3.1f}.txt".format(temperature_k[k]))
    U_kn[k, :] = lines[1500:]

U_kn = np.array(U_kn)

avgE = np.mean(U_kn, axis=1)
# print avgE.shape

plt.figure(0)
plt.plot(temperature_k, avgE, 'o')
plt.show()
