import numpy as np
import matplotlib.pyplot as plt

def read_file(filename):
    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()
    return lines

temperature_k = [349.0, 349.5, 350.0, 350.5, 351.0, 351.5]
K = len(temperature_k)
T = 1001

E_kn = np.zeros((K, T))		# U_kn[k,t] is the potential energy (in kcal/mol) for snapshot t of temperature index k
for k in range(K):
    lines = read_file("en-new.{:3.1f}.txt".format(temperature_k[k]))
    E_kn[k, :] = lines[1000:]

E_kn = np.array(E_kn)

avgE = np.mean(E_kn, axis=1)
# print avgE.shape

plt.figure(0)
plt.plot(temperature_k, avgE, 'o', markersize=20)
plt.show()


theta_kn = np.zeros((K, T))		# U_kn[k,t] is the potential energy (in kcal/mol) for snapshot t of temperature index k
for k in range(K):
    theta_k = np.genfromtxt("theta{:3.1f}.txt".format(temperature_k[k]))
    theta_k = theta_k.reshape((-1, 240))
    theta_lines = np.mean(theta_k, axis=1)
    theta_kn[k, :] = theta_lines[1000:]

theta_kn = np.array(theta_kn)

avgtheta = np.mean(theta_kn, axis=1)
# print avgtheta.shape

plt.figure(1)
plt.plot(temperature_k, avgtheta, 'o', markersize=20)
plt.show()
