import numpy as np
import pandas as pd
from scipy.integrate import odeint
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from matplotlib import style

style.use('default')

data = pd.read_csv(
    '/Users/AymenHafeez/Desktop/DESKTOP/Research_project/Microalgae_Models/data_4.15.csv')

t1 = data['Day'].values
X1meas = data['dry_weight2'].values
ns = len(t1)


def model(x, t, p):
    X1 = x[0]
    sl = x[1]
    sg = x[2]

    muMax, ks, ki, k, Yi, i0, a = p
    i = i0 / (a * X1 * (1 - np.exp(-a * X1)))
    mu = muMax * sl / (sl + ks + (sl**2) / ki * (i / (i + k)))
    kla = 0.00095
    h = 0.00316
    sgin = 0.06

    dXdt = mu * X1
    dsldt = kla * ((sg / h) - sl) - (Yi * X1)
    dsgdt = sgin - kla * ((sg / h) - sl)
    return [dXdt, dsldt, dsgdt]


X0 = 0.028
sl0 = 0.0
sg0 = 0.0


def simulate(p):
    X = np.zeros((len(t1), 3))
    X[0] = X1meas[0]
    X0 = X[0]
    for i in range(len(t1) - 1):
        ts = [t1[i], t1[i + 1]]
        x = odeint(model, X0, ts, args=(p, ))
        X0 = x[-1]
        X[i + 1] = X0
    return X


def objective(p):
    Xp = simulate(p)
    obj = 0.0
    for i in range(len(t1)):
        obj += ((Xp[i, 0] - X1meas[i]) / X1meas[i])**2
    return obj


# Initial parameter values
muMax = 0.8
ks = 0.03
ki = 3
k = 14
Yi = 0.05
i0 = 75
a = 0.014
p0 = [muMax, ks, ki, k, Yi, i0, a]

solution = minimize(objective, p0, method='SLSQP')
p = solution.x

print('Initial SSE: ' + str(objective(p0)))
print('Final SSE: ' + str(objective(p)))

# Optimised parameter values
muMax = p[0]
ks = p[1]
ki = p[2]
k = p[3]
Yi = p[4]
i0 = p[5]
a = p[6]

print('muMax: ' + str(abs(muMax)))
print('K_s: ' + str(ks))
print('K_i: ' + str(ki))
print('K: ' + str(k))
print('Y_i: ' + str(Yi))
print('i_0: ' + str(i0))
print('a: ' + str(a))

xi = simulate(p0)
xp = simulate(p)

plt.plot(t1, X1meas, '.', markersize=8, markerfacecolor='w',
         markeredgecolor='b', markeredgewidth=0.5, label='Measured data')
plt.plot(t1, xi[:, 0], 'g--', linewidth=1, label='Initial prediction')
plt.plot(t1, xi[:, 1], 'r-', linewidth=1, label='Optimised model')
plt.ylabel(r'Biomass concentration (g L$^{-1}$)')
plt.xlabel('Time (days)')
plt.legend(fontsize=12)
plt.ylim(-0.2, 2.75)

plt.show()
