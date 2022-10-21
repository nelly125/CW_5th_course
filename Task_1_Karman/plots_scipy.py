import numpy as np
import matplotlib.pyplot as plt
import pylab
import scipy.integrate

s = 0.5

def f(ksi, y):
    F = y[0]
    R = y[1]
    G = y[2]
    T = y[3]
    H = y[4]

    dF = R
    dR = H * R + F * F - G * G + s * s
    dG = T
    dT = 2 * F * G + H * T
    dH = -2 * F
    return [dF, dR, dG, dT, dH]


initial_condition = [0,  0.391, 1, -0.4176, 0]
# print(f(0, [1, 1, 25, 1, 1]));
left_boarder = 0
right_boarder = 5

res = scipy.integrate.solve_ivp(fun=f, t_span=[0., 75.], y0=initial_condition, method='RK45',
                                t_eval=np.linspace(left_boarder, right_boarder, 1000))

fig = plt.figure(figsize=(20, 15))
plt.plot(res.t, res.y[0], linewidth=4, markersize=1.0, label="F")
plt.plot(res.t, res.y[2], linewidth=4, markersize=1.0, label="G")
plt.plot(res.t, res.y[4], linewidth=4, markersize=1.0, label="H")
plt.title("Parameters: " + "alpha: " + str(initial_condition[1]) + " beta: " + str(initial_condition[3]), fontsize=20)

plt.xlabel("ksi", fontsize=20)
plt.ylabel("Solution", fontsize=20)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(fontsize=16)
plt.grid()
plt.show()
plt.savefig("../RK_4_scipy_plots.png")
