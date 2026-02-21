"""
MAE 6760 Model Based Estimation
Cornell University
M Campbell

Homework #2
Problem #5: Nonlinear Least Squares
for multi-sensor range localization problem
"""

import numpy as np

np.set_printoptions(precision=5, suppress=True)
np.random.seed(10)

# =================================================
# Problem setup: beacons and measurements
# =================================================
bA = np.array([[-10.0], [100.0]])
bB = np.array([[490.0], [20.0]])
bC = np.array([[500.0], [40.0]])
beacons = np.hstack((bA, bB, bC))

xtrue = np.array([[-5.0], [2.0]])

# Perfect ranges
RA = np.linalg.norm(bA - xtrue)
RB = np.linalg.norm(bB - xtrue)
RC = np.linalg.norm(bC - xtrue)

n = 10  # number of measurements

# =================================================
# Part (a): three beacons (A,B,C)
# =================================================
Rpart_a = np.diag([10.0, 10.0, 10.0])
v_a = np.linalg.cholesky(Rpart_a) @ np.random.randn(3, n)
z_a = np.tile(np.array([[RA], [RB], [RC]]), (1, n)) + v_a


# =================================================
# Part (b): two beacons (A,B) and (B,C)
# =================================================
# (A,B)
ii = [0, 1]
Rpart_b1 = Rpart_a[np.ix_(ii, ii)]
z_b1 = z_a[ii, :]


# (B,C)
ii = [1, 2]
Rpart_b2 = Rpart_a[np.ix_(ii, ii)]
z_b2 = z_a[ii, :]


# =================================================
# Part (c): perfect linearization
# =================================================


# =================================================
# Part (d): correlated noise
# =================================================
Rpart_d = np.array([[10, 0, 0],
                    [0, 10, 9],
                    [0, 9, 10]], dtype=float)

v_d = np.linalg.cholesky(Rpart_d) @ np.random.randn(3, n)
z_d = np.tile(np.array([[RA], [RB], [RC]]), (1, n)) + v_d

