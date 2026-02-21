"""
MAE 6760 Model Based Estimation
Cornell University
M Campbell

Homework #2
Problem #4: MMSE Estimation
for a time-varying output function
"""

import numpy as np
import matplotlib.pyplot as plt

np.set_printoptions(precision=8, suppress=True)

# =================================================
# Simulate the measurements
# =================================================
t = np.arange(0.1, 1.0 + 0.1, 0.1).reshape(-1, 1)
n = len(t)

x = np.array([[1.0],
              [1.0],
              [1.0]])

R_diag = np.array([0.001, 0.002, 0.005, 0.010, 0.008,
                   0.002, 0.010, 0.007, 0.020, 0.006])**2
R = np.diag(R_diag)

v = np.linalg.cholesky(R) @ np.random.randn(n, 1)

z = (
    x[0]
    + x[1] * np.sin(10 * t)
    + x[2] * np.exp(2 * t**2)
    + v
)

# =================================================
# Prior information
# =================================================
x0 = np.array([[1.01],
               [0.98],
               [0.99]])

P0 = np.eye(3) * 0.001

# =================================================
# Part (a): MMSE estimate
# =================================================


# =================================================
# Part (b): Varying sensor noise
# =================================================

# plotting routine for error bars, assuming 3x7 xi and sigi

# Plot results
fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharex=True)
true_val = np.ones(ni)

for i in range(3):
    axes[i].errorbar(alphai, xi[i, :], 2 * sigi[i, :], fmt='b')
    axes[i].semilogx(alphai, true_val, color=(200/255, 0, 0))
    axes[i].set_xscale('log')
    axes[i].grid(True)
    axes[i].set_xlabel('scale factor α')
    axes[i].set_ylabel('state estimate')

axes[1].legend(['est ± 2σ', 'true value'], loc='lower left')
fig.suptitle(r'study of sensor noise scale factor $\alpha$ ($R=\alpha R_0$)')

# =================================================
# Part (c): Varying prior
# =================================================

# plotting routine for error bars, assuming 3x7 xi and sigi

# Plot results
fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharex=True)

for i in range(3):
    axes[i].errorbar(alphai, xi[i, :], 2 * sigi[i, :], fmt='b')
    axes[i].semilogx(alphai, true_val, color=(200/255, 0, 0))
    axes[i].set_xscale('log')
    axes[i].grid(True)
    axes[i].set_xlabel('scale factor α')
    axes[i].set_ylabel('state estimate')

axes[1].legend(['est ± 2σ', 'true value'], loc='lower left')
fig.suptitle(r'study of prior scale factor $\alpha$ ($P=\alpha P_0$)')

plt.show()