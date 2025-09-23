import numpy as np
import matplotlib.pyplot as plt
# Chatgpt was used in helping write this code
# --------------------------
# Helper functions
# --------------------------
def ranges(p, B):
    return np.linalg.norm(B - p[None, :], axis=1)

def jacobian_range_wrt_p(p, B):
    d = B - p[None, :]
    r = np.linalg.norm(d, axis=1, keepdims=True)
    r = np.maximum(r, 1e-9)
    return (B - p) / r

def pdop_from_geometry(p, B):
    J = jacobian_range_wrt_p(p, B)
    JTJ = J.T @ J
    return np.sqrt(np.trace(np.linalg.inv(JTJ)))

def gn_solve_position(B, r_meas, p0, iters=15):
    p = p0.copy()
    for _ in range(iters):
        r = ranges(p, B)
        e = r - r_meas
        u = (p[None, :] - B) / np.maximum(r[:, None], 1e-9)
        JTJ = u.T @ u
        JTe = u.T @ e
        dp = -np.linalg.solve(JTJ, JTe)
        p = p + dp
        if np.linalg.norm(dp) < 1e-9:
            break
    return p

def map_solve_with_beacon_uncertainty(B_mu, r_meas, sigma_b=0.0, p0=None, iters=20):
    p = p0.copy()
    B = B_mu.copy()
    sig2_b = sigma_b**2
    for _ in range(iters):
        # Step A: update p
        p = gn_solve_position(B, r_meas, p, iters=5)
        # Step B: update B with priors
        if sig2_b > 0:
            for i in range(B.shape[0]):
                Bi = B[i].copy()
                mu = B_mu[i]
                di = Bi - p
                ri = np.linalg.norm(di)
                ri = max(ri, 1e-9)
                ui = di / ri
                ei = ri - r_meas[i]
                A = np.outer(ui, ui) + (1.0/sig2_b) * np.eye(2)
                b = -(ei * ui + (1.0/sig2_b) * (Bi - mu))
                dBi = np.linalg.solve(A, b)
                B[i] = Bi + dBi
    return p, B

# --------------------------
# Setup
# --------------------------
p_true = np.array([8.846, 3.740])  # true SCOUT position at t=10s
B_mu = np.array([
    [ 100.0, -100.0],
    [-100.0,  100.0],
    [ 100.0,  100.0],
    [-100.0, -100.0],
])

rng = np.random.default_rng(7)
r_true = ranges(p_true, B_mu)
multipliers = rng.choice([0.95, 1.05], size=4)
r_meas = r_true * multipliers

# Case A: fixed beacons
p0 = B_mu.mean(axis=0)
p_est_fixed = gn_solve_position(B_mu, r_meas, p0)

# Case B: beacon uncertainty
sigma_b = 2.0
p_est_map, B_est = map_solve_with_beacon_uncertainty(B_mu, r_meas, sigma_b, p0)

# --------------------------
# PDOP values
# --------------------------
pdop_truth = pdop_from_geometry(p_true, B_mu)
pdop_fixed = pdop_from_geometry(p_est_fixed, B_mu)
pdop_uncertain = pdop_from_geometry(p_est_map, B_est)

print("PDOP at true position:       ", pdop_truth)
print("PDOP at est. (fixed beacons):", pdop_fixed)
print("PDOP at est. (σ_b=2m):       ", pdop_uncertain)
sigmas_b = [0.0, 1.0, 2.0, 5.0]  # m
rows = []
for sb in sigmas_b:
    if sb == 0.0:
        # fixed beacons: estimate p with nominal beacon positions
        p_est = gn_solve_position(B_mu, r_meas, p0)          # your GN solver
        B_use = B_mu
    else:
        # draw shifted beacons from N(B_mu, sb^2 I) and re-solve
        B_shift = B_mu + rng.normal(0.0, sb, size=B_mu.shape)
        p_est   = gn_solve_position(B_shift, r_meas, p0)     # or your MAP joint solver
        B_use   = B_shift

    # PDOP at the estimate (geometry only)
    pd = pdop_from_geometry(p_est, B_use)

    # error vs truth
    err = np.linalg.norm(p_est - p_true)

    rows.append((sb, p_est[0], p_est[1], err, pd))

print("sigma_b[m]    x_hat[m]    y_hat[m]   pos_err[m]   PDOP")
for sb,xh,yh,err,pd in rows:
    print(f"{sb:6.2f}    {xh:8.2f}  {yh:8.2f}    {err:8.2f}   {pd:.6f}")

# --------------------------
# Plot
# --------------------------
plt.figure(figsize=(7,6))
plt.scatter(B_mu[:,0], B_mu[:,1], c='k', marker='s', s=80, label="Nominal beacons")
plt.scatter(B_est[:,0], B_est[:,1], c='orange', marker='s', s=80, label=f"Shifted beacons (σ_b={sigma_b} m)")
plt.scatter(*p_true, c='g', s=80, label="True AUV pos")
plt.scatter(*p_est_fixed, c='b', s=80, label="Est. pos (fixed beacons)")
plt.scatter(*p_est_map, c='r', s=80, label=f"Est. pos (σ_b={sigma_b} m)")
for mu, est in zip(B_mu, B_est):
    plt.plot([mu[0], est[0]], [mu[1], est[1]], 'k--', alpha=0.5)
plt.axis('equal')
plt.grid(True)
plt.xlabel("X [m]")
plt.ylabel("Y [m]")
plt.title("LBL Position Estimation with Beacon Uncertainty")
plt.legend()
plt.tight_layout()
plt.show()
