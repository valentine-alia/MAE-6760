"""
Sigma Point Filter Example
    planar pose of a car (four states)
    nonlinear dynamics, IMU and 2D position (GPS-like) measurements

Requires:
    plot_estimator.py
    calculate_ellipse.py

MAE 6760 Model Based Estimation
Cornell University
M Campbell
"""

import numpy as np
import matplotlib.pyplot as plt

from plot_estimator import plot_estimator
from calculate_ellipse import calculate_ellipse

# ----------------------------------------------------------
# USER INPUT
# ----------------------------------------------------------

nsig = 1

# ----------------------------------------------------------
# DEFINE SYSTEM MODEL AND SIMULATE
# ----------------------------------------------------------

np.random.seed(0)

Uacc = np.concatenate(
    [
        np.ones(60),
        np.zeros(100),
        -np.ones(60),
        np.zeros(20),
        np.ones(60),
        np.zeros(50),
    ]
)

Uomegadot = np.concatenate(
    [
        np.zeros(240),
        -np.ones(40) * np.pi / 2 / 4,
        np.zeros(70),
    ]
)

dt = 0.1
nt = len(Uacc)
t = np.arange(0, nt * dt, dt)

nx = 4
x_true = np.zeros((nx, nt))

for k in range(nt - 1):

    Vk = x_true[2, k]
    Tk = x_true[3, k]

    x_true[:, k + 1] = x_true[:, k] + dt * np.array(
        [
            Vk * np.cos(Tk),
            Vk * np.sin(Tk),
            Uacc[k],
            Uomegadot[k],
        ]
    )

# ----------------------------------------------------------
# PLOT TRUE TRAJECTORY
# ----------------------------------------------------------

plt.figure(figsize=(8,6))
plt.plot(
    x_true[0, :],
    x_true[1, :],
    "m.",
)
plt.axis([0, 110, -50, 10])
plt.grid()
plt.ylabel("North (m)")
plt.xlabel("East (m)")
plt.tight_layout()

# ----------------------------------------------------------

plt.figure(figsize=(8,6))

plt.plot(t, x_true[2, :], "b-")
plt.plot(t, x_true[3, :], "b--")

plt.ylabel("ideal state")
plt.xlabel("time (sec)")
plt.legend(["velocity (m/sec)", "heading (rad)"], loc="lower left")
plt.axis([0, nt * dt, -2, 6.5])
plt.grid()
plt.tight_layout()

# ----------------------------------------------------------
# PROCESS NOISE AND IMU MEASUREMENTS
# ----------------------------------------------------------

nw = 2
Q = np.diag([(dt * np.max(np.abs(Uacc))) ** 2, (dt * np.max(np.abs(Uomegadot))) ** 2])
w = np.linalg.cholesky(Q) @ np.random.randn(nw, nt)

Zacc = Uacc + w[0, :]
Zomegadot = Uomegadot + w[1, :]

G = np.array(
    [
        [0, 0],
        [0, 0],
        [dt, 0],
        [0, dt],
    ]
)
Qx = G @ Q @ G.T

# ----------------------------------------------------------
# GPS MEASUREMENTS
# ----------------------------------------------------------

nz = 2
R = np.eye(nz) * 1**2
v = np.linalg.cholesky(R) @ np.random.randn(nz, nt)
z = x_true[0:2, :] + v
H = np.hstack((np.eye(2), np.zeros((2, 2))))

# ----------------------------------------------------------
# FUNCTIONS
# ----------------------------------------------------------

def predict_state_carpose(Xk, U, dt):

    Vk = Xk[2, :]
    Tk = Xk[3, :]

    acc = U[0] * np.ones_like(Vk)
    omegadot = U[1] * np.ones_like(Vk)

    Xkp1 = Xk + dt * np.vstack(
        [
            Vk * np.cos(Tk),
            Vk * np.sin(Tk),
            acc,
            omegadot,
        ]
    )

    return Xkp1

def predict_msmt_carpose(Xkp1):

    H = np.array(
        [
            [1, 0, 0, 0],
            [0, 1, 0, 0],
        ]
    )

    return H @ Xkp1

def spfadd(xEst, PxEst, U, Q, z, R, dt, nsig):

    nx = len(xEst)
    nsp = 2 * nx + 1

    Wi = 0.5 / nsig**2
    W0M = (nsig**2 - nx) / nsig**2
    W0C = (nsig**2 - nx) / nsig**2 + 3 - nsig**2 / nx

    WM = np.concatenate(([W0M], np.ones(2 * nx) * Wi))
    WC = np.concatenate(([W0C], np.ones(2 * nx) * Wi))
    # print(PxEst)
    Psqrtm = nsig * np.linalg.cholesky(PxEst)

    xSigmaPts = np.hstack(
        (
            np.zeros((nx, 1)),
            -Psqrtm,
            Psqrtm,
        )
    )

    xSigmaPts = xSigmaPts + xEst[:, None]
    xPredSigmaPts = predict_state_carpose(xSigmaPts, U, dt)
    xPred = (xPredSigmaPts @ WM).flatten()

    exSigmaPts = xPredSigmaPts - xPred[:, None]
    PxxPred = (exSigmaPts * WC) @ exSigmaPts.T + Q

    Psqrtm = nsig * np.linalg.cholesky(PxxPred)
    exSigmaPts = np.hstack((np.zeros((nx, 1)), -Psqrtm, Psqrtm))
    xPredSigmaPts = exSigmaPts + xPred[:, None]
    zPredSigmaPts = predict_msmt_carpose(xPredSigmaPts)

    zPred = (zPredSigmaPts @ WM).flatten()
    ezSigmaPts = zPredSigmaPts - zPred[:, None]

    PxzPred = (exSigmaPts * WC) @ ezSigmaPts.T
    Pyy = (ezSigmaPts * WC) @ ezSigmaPts.T + R

    K = PxzPred @ np.linalg.inv(Pyy)
    PxEst = PxxPred - K @ PxzPred.T
    # PxEst = PxxPred - K @ Pyy @ K.T
    PxEst = 0.5 * (PxEst + PxEst.T)

    innovation = z - zPred
    xEst = xPred + K @ innovation

    return xEst, PxEst, xPred, zPred, innovation


# ----------------------------------------------------------
# SPF INITIALIZATION
# ----------------------------------------------------------

x0 = np.array([0, 0, 0, 0])

P0 = np.diag([2**2, 2**2, 1**2, 0.1**2])

xhatu = np.zeros((nx, nt))
xhatu[:, 0] = x0

Pu = np.zeros((nx, nx, nt))
Pu[:, :, 0] = P0

# ----------------------------------------------------------
# SPF LOOP
# ----------------------------------------------------------

for k in range(nt - 1):

    Uk = np.array([Zacc[k], Zomegadot[k]])

    (
        xhatu[:, k + 1],
        Pu[:, :, k + 1],
        xPred,
        zPred,
        innovation,
    ) = spfadd(
        xhatu[:, k],
        Pu[:, :, k],
        Uk,
        Qx,
        z[:, k + 1],
        R,
        dt,
        nsig,
    )

# ----------------------------------------------------------
# PLOTS
# ----------------------------------------------------------

fig = plt.figure(figsize=(16, 6))
gs = fig.add_gridspec(1, 2)

ax = fig.add_subplot(gs[0])
plot_estimator(
    t,
    xhatu[0, :],
    Pu[0, 0, :],
    x_true[0, :],
    "error",
    z[0, :],
    ax=ax,
)
ax.set_ylabel("x position (m)")
ax.axis([0, 35, -3, 3])

ax = fig.add_subplot(gs[1])
plot_estimator(
    t,
    xhatu[1, :],
    Pu[1, 1, :],
    x_true[1, :],
    "error",
    z[1, :],
    ax=ax,
)
ax.set_ylabel("y position (m)")
ax.axis([0, 35, -3, 3])
plt.tight_layout()


# ----------------------------------------------------------

fig = plt.figure(figsize=(16, 6))
gs = fig.add_gridspec(1, 2)

ax = fig.add_subplot(gs[0])
plot_estimator(
    t,
    xhatu[2, :],
    Pu[2, 2, :],
    x_true[2, :],
    "error",
    ax=ax,
)
ax.set_ylabel("V velocity (m/sec)")
ax.axis([0, 35, -1, 1])

ax = fig.add_subplot(gs[1])

thC = 180 / np.pi

plot_estimator(
    t,
    xhatu[3, :] * thC,
    Pu[3, 3, :] * thC**2,
    x_true[3, :] * thC,
    "error",
    ax=ax,
)

ax.set_ylabel(r"$\theta$ heading (deg)")
ax.axis([0, 35, -20, 20])
plt.tight_layout()

# ---------------------------------------------------------
# Birds eye view trajectory with covariance ellipses
# ---------------------------------------------------------

plt.figure(figsize=(8,6))

pt, = plt.plot(x_true[0,:],x_true[1,:],'m.',linewidth=3)
ph, = plt.plot(xhatu[0,:],xhatu[1,:],'b-')
iell = np.concatenate(([2],np.arange(10,nt,10)))
patch_handles=[]

for ii in iell:
    Xe,Ye,_,_,_ = calculate_ellipse(
        xhatu[0:2,ii],
        Pu[0:2,0:2,ii],
        3,
        50
    )
    plt.plot(Xe, Ye, 'b-', linewidth=0.5)
    plt.plot(xhatu[0,ii],xhatu[1,ii],'bx')
    p = plt.fill(Xe,Ye,'b',alpha=0.1)
    patch_handles.append(p[0])

plt.xlabel('x position (m)')
plt.ylabel('y position (m)')
plt.legend([ph,patch_handles[0],pt],
           ['position estimate','error ellipse','truth'],
           loc='lower left')
plt.axis([-10,110,-50,10])
plt.tight_layout()

plt.show()