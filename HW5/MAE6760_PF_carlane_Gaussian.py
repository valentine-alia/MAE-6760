"""
Particle Filter Example:
  planar pose of a car (four states)
  nonlinear dynamics
  odometry and lane based position measurements

  Note: uses plot_estimator.py, calculate_ellipse.py

MAE 6760 Model Based Estimation
Cornell University
M Campbell
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from plot_estimator import plot_estimator
from calculate_ellipse import calculate_ellipse

# =============================================================================
# User parameters
# =============================================================================
N_part      = 1000
Neff_thresh = N_part / 2   # common approach
# Neff_thresh = 1           # no resampling: SIS
# Neff_thresh = N_part      # always resample: SIR

# =============================================================================
# Reproducibility
# =============================================================================
rng = np.random.default_rng(0)

# =============================================================================
# MAE6760startup — figure / axes defaults
# =============================================================================
# plt.rcParams.update({
#     "figure.figsize":   (8, 6),
#     "figure.dpi":       100,
#     "axes.titlesize":   14,
#     "axes.labelsize":   14,
#     "axes.titleweight": "bold",
#     "axes.labelweight": "bold",
#     "xtick.labelsize":  12,
#     "ytick.labelsize":  12,
#     "lines.linewidth":  2,
#     "font.weight":      "bold",
# })

# =============================================================================
# Internal functions
# =============================================================================

def predict_state_carpose(Xk, Wk, dt):
    """
    For the car localization/pose example.
    Discrete prediction of state at k+1.

    Constant velocity/heading model assumed (no external inputs).

    Parameters
    ----------
    Xk  : (4, N) array  — state [x; y; V; Theta]
    Wk  : (2, N) array  — process noise [W_acc; W_rg]
    dt  : float          — time step

    Returns
    -------
    Xkp1 : (4, N) array
    """
    Vk    = Xk[2, :]
    Tk    = Xk[3, :]
    Wacck = Wk[0, :]
    Wrgk  = Wk[1, :]

    Xkp1 = Xk + dt * np.vstack([
        Vk * np.cos(Tk),
        Vk * np.sin(Tk),
        Wacck,
        Wrgk,
    ])
    return Xkp1

# =============================================================================
# Set-up model, simulate true car position
# =============================================================================
nx = 4          # state X = [x, y, V, Theta]
lw = 3.7        # lane width (m)
x0 = np.array([-lw / 2, 0.0, 5.0, np.pi / 2])

dt = 0.1
nk = 100
t  = np.arange(nk) * dt    # (nk,)

x_true = np.zeros((nx, nk))
x_true[:, 0] = x0
for k in range(nk - 1):
    x_true[:, k + 1] = predict_state_carpose(
        x_true[:, k:k + 1], np.zeros((2, 1)), dt
    ).ravel()

x = x_true[0, :]
y = x_true[1, :]

# =============================================================================
# Sensor model: lane detector
# =============================================================================
vx_sig = 0.5
vx_mn  = 0.0
vx     = rng.normal(vx_mn, vx_sig, nk)   # lane x measurement noise

# =============================================================================
# Sensor model: odometry
# =============================================================================
vy_sig = 0.5
vy_mn  = 0.0
vy     = rng.normal(vy_mn, vy_sig, nk)   # odometry y measurement noise

# =============================================================================
# Full measurement z and process noise matrices
# =============================================================================
nz = 2
H  = np.hstack([np.eye(2), np.zeros((2, 2))])
z  = H @ x_true + np.vstack([vx, vy])    # (2, nk)
R  = np.diag([vx_sig**2, vy_sig**2])

nw  = 2
Q   = np.diag([0.5**2, 0.1**2])
Qsq = np.linalg.cholesky(Q)              # square root of Q (lower triangular)

# =============================================================================
# Initialise state particles and weights
# =============================================================================
P0     = np.diag([0.5**2, 0.5**2, 0.5**2, 0.1**2])
X_part = np.zeros((nx, N_part, nk))
X_part[:, :, 0] = rng.multivariate_normal(x0, P0, N_part).T

W0 = np.ones(N_part) / N_part            # uniform initial weights
W  = np.zeros((nk, N_part))
W[0, :] = W0

# =============================================================================
# Initial weighted-sample estimate and covariance
# =============================================================================
xhat = np.zeros((nx, nk))
xhat[:, 0] = (W0 * X_part[:, :, 0]).sum(axis=1)

Err0 = X_part[:, :, 0] - xhat[:, 0:1]
Phat = np.zeros((nx, nx, nk))
Phat[:, :, 0] = (Err0 * W0) @ Err0.T

# =============================================================================
# Figure 10: lane view + zoomed particle panel
# =============================================================================
ylen  = 52
fig10 = plt.figure(10, figsize=(10, 8))
gs    = GridSpec(4, 4, figure=fig10)

# Left panel  — mirrors subplot(4,4,[1 5 9 13])
ax_left  = fig10.add_subplot(gs[:, 0:2])
# Right panel — mirrors subplot(4,4,[7:8 11:12])
ax_right = fig10.add_subplot(gs[1:3, 2:4])

for ax in (ax_left, ax_right):
    ax.axvline(0,       color="k",    linewidth=3)
    ax.axvline(-2 * lw, color="k",    linewidth=3)
    ax.axvline(-lw,     color="gold", linewidth=5, linestyle="--")
    ax.set_facecolor([0.5, 0.5, 0.5])
    ax.set_xlabel("cross-lane $x$")
    ax.set_ylabel("along-lane $y$")
    ax.set_title("particle distribution")

ax_left.set_xticks([-2 * lw, -lw, 0])
ax_left.set_xlim(-2 * lw, 0)
ax_left.set_ylim(-2, ylen)

ax_right.set_xticks([-2 * lw, -1.5 * lw, -lw, -0.5 * lw, 0])
ax_right.set_xlim(-7.4, 0)
ax_right.set_ylim(-4, 4)

# Initial ellipse
Xe0, Ye0, _, _, _ = calculate_ellipse(xhat[0:2, 0], Phat[0:2, 0:2, 0], nsig=2)

# Left panel handles
h_pts,  = ax_left.plot(X_part[0, :, 0], X_part[1, :, 0], "b.", markersize=6)
h_est,  = ax_left.plot(xhat[0, 0], xhat[1, 0], "m*", markersize=10, linewidth=2)
h_ell,  = ax_left.plot(Xe0, Ye0, "m-", linewidth=2)

# Right panel handles
h_pts2, = ax_right.plot(X_part[0, :, 0], X_part[1, :, 0], "b.", markersize=6)
h_est2, = ax_right.plot(xhat[0, 0], xhat[1, 0], "m*", markersize=10, linewidth=2)
h_ell2, = ax_right.plot(Xe0, Ye0, "m-", linewidth=2)

plt.tight_layout()
plt.ion()
plt.show()

# =============================================================================
# Run particle filter
# =============================================================================
Neff    = np.zeros(nk)
Neff[0] = N_part
Wk      = W0.copy()
CDF     = np.cumsum(Wk)
L       = np.zeros(N_part)
R_inv   = np.linalg.inv(R)

for k in range(nk - 1):

    X_prior = X_part[:, :, k]

    # ------------------------------------------------------------------
    # PREDICTION STEP
    # ------------------------------------------------------------------
    W_part = Qsq @ rng.standard_normal((nw, N_part))   # sample process noise
    X_pred = predict_state_carpose(X_prior, W_part, dt)

    # ------------------------------------------------------------------
    # UPDATE STEP
    # ------------------------------------------------------------------
    z_current = z[:, k + 1]
    Z_hat     = H @ X_pred                              # estimated measurements (2, N)
    Inn       = z_current[:, None] - Z_hat              # innovations (2, N)

    # Likelihood weighting for each particle
    for ip in range(N_part):
        inn_ip = Inn[:, ip]
        L[ip]  = np.exp(-0.5 * inn_ip @ R_inv @ inn_ip)

    # Update weights
    Wk_unnorm = W[k, :] * L
    Wk        = Wk_unnorm / Wk_unnorm.sum()
    CDF       = np.cumsum(Wk) / Wk.sum()

    # ------------------------------------------------------------------
    # RESAMPLING
    # ------------------------------------------------------------------
    Neff[k] = 1.0 / (Wk @ Wk)    # effective number of particles

    if Neff[k] < Neff_thresh:
        CDF      = np.cumsum(Wk) / Wk.sum()
        CDF_plus = CDF + rng.uniform(0, 1e-6, N_part)  # jitter for zero-weight particles
        iSelect  = rng.uniform(0, 1, N_part)
        iNextGen = np.clip(np.searchsorted(CDF_plus, iSelect), 0, N_part - 1)
        X_part[:, :, k + 1] = X_pred[:, iNextGen]
        W[k + 1, :]          = np.ones(N_part) / N_part
        Wk_plot              = np.ones(N_part) / N_part
    else:
        X_part[:, :, k + 1] = X_pred
        W[k + 1, :]          = Wk
        Wk_plot              = Wk

    # Weighted mean estimate and sample covariance
    xhat[:, k + 1]    = (Wk_plot * X_part[:, :, k + 1]).sum(axis=1)
    Err                = X_part[:, :, k + 1] - xhat[:, k + 1:k + 2]
    Phat[:, :, k + 1] = (Err * Wk_plot) @ Err.T

    # ------------------------------------------------------------------
    # Redraw particles, estimate, ellipse
    # ------------------------------------------------------------------
    Xe, Ye, _, _, _ = calculate_ellipse(xhat[0:2, k + 1], Phat[0:2, 0:2, k + 1], nsig=2)

    h_pts.set_xdata(X_part[0, :, k + 1]);  h_pts.set_ydata(X_part[1, :, k + 1])
    h_est.set_xdata([xhat[0, k + 1]]);     h_est.set_ydata([xhat[1, k + 1]])
    h_ell.set_xdata(Xe);               h_ell.set_ydata(Ye)

    h_pts2.set_xdata(X_part[0, :, k + 1]); h_pts2.set_ydata(X_part[1, :, k + 1])
    h_est2.set_xdata([xhat[0, k + 1]]);    h_est2.set_ydata([xhat[1, k + 1]])
    h_ell2.set_xdata(Xe);              h_ell2.set_ydata(Ye)
    ax_right.set_ylim(y[k + 1] - 4, y[k + 1] + 4)

    fig10.canvas.draw_idle()
    plt.pause(0.1)

plt.ioff()

# =============================================================================
# Extra plots
# =============================================================================

# --- state estimate + 2-sigma bounds ---
fig1, axes1 = plt.subplots(2, 1, figsize=(8, 6))
plot_estimator(t, xhat[0, :], Phat[0, 0, :], x_true[0, :],
               plot_type="state", z=z[0, :], ax=axes1[0])
axes1[0].set_ylabel("x position (m)")
axes1[0].set_yticks([-2 * lw, -lw, 0])
axes1[0].set_ylim(-2 * lw, 0)

plot_estimator(t, xhat[1, :], Phat[1, 1, :], x_true[1, :],
               plot_type="state", z=z[1, :], ax=axes1[1])
axes1[1].set_ylabel("y position (m)")
plt.tight_layout()

# --- error + 2-sigma bounds ---
fig2, axes2 = plt.subplots(2, 1, figsize=(8, 6))
plot_estimator(t, xhat[0, :], Phat[0, 0, :], x_true[0, :],
               plot_type="error", z=z[0, :], ax=axes2[0])
axes2[0].set_ylim(-1.5, 1.5)
axes2[0].set_ylabel("error x position (m)")
axes2[0].get_legend().remove()

plot_estimator(t, xhat[1, :], Phat[1, 1, :], x_true[1, :],
               plot_type="error", z=z[1, :], ax=axes2[1])
axes2[1].set_ylim(-1.5, 1.5)
axes2[1].set_ylabel("error y position (m)")
plt.tight_layout()

# --- number of effective particles ---
fig3, ax3 = plt.subplots(figsize=(8, 6))
ax3.semilogy(np.arange(1, nk + 1), Neff,
             "-", label="# of effective particles")
ax3.semilogy(np.arange(1, nk + 1), Neff_thresh * np.ones(nk),
             "r:", linewidth=3, label="threshold")
ax3.set_ylabel("# of effective particles")
ax3.set_xlabel("timestep k")
ax3.legend()
ax3.set_xlim(0, 100)
ax3.set_ylim(1, 1000)
ax3.grid(True)
plt.tight_layout()

# --- final particle weights ---
fig4, ax4 = plt.subplots(figsize=(8, 6))
ax4.semilogy(np.arange(1, N_part + 1), Wk, "mx")
ax4.set_ylabel("particle weights")
ax4.set_xlabel("particle number")
ax4.set_xlim(0, N_part)
ax4.set_ylim(0, 0.1)
plt.tight_layout()

# --- CDF of particle weights ---
fig5, ax5 = plt.subplots(figsize=(8, 6))
ax5.plot(np.arange(1, N_part + 1), CDF, "-")
ax5.set_ylabel("CDF of particle weights")
ax5.set_xlabel("particle number")
plt.tight_layout()

plt.show()