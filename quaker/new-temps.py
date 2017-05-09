import numpy as np
import matplotlib.pyplot as plt

log = np.genfromtxt('log.stripped').astype(int)

T_n_up = np.zeros(6)
T_n_down = np.zeros(6)
R_uplabels = [None, None, None, None, None, None]

for T_from_R in log:
    zero_traj = np.where( T_from_R[1:] == 0 )[0][0]
    N_traj = np.where( T_from_R[1:] == 5 )[0][0]
#     print zero_traj, N_traj
    R_uplabels[zero_traj] = True
    R_uplabels[N_traj] = False
    for R_i in range(6):
	T_i = T_from_R[R_i+1]
        if (R_uplabels[R_i] is None):
            pass
        elif (R_uplabels[R_i]):
            T_n_up[T_i] = T_n_up[T_i] + 1
        else:
            T_n_down[T_i] = T_n_down[T_i] + 1

print T_n_up, T_n_down

T_f_up = T_n_up / (T_n_up + T_n_down)

plt.figure(0)
plt.plot(T_f_up)
plt.show()

