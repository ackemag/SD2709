#!/usr/bin/env python3
# Kalman Filter for SCOUT AUV (constant-velocity model with LBL updates)
# Chatgpt was used in helping write this code
import numpy as np
import matplotlib.pyplot as plt

def build_cv_matrices(dt, sigma_a, sigma_meas):
    F = np.array([[1, 0, dt, 0],
                  [0, 1, 0, dt],
                  [0, 0, 1, 0],
                  [0, 0, 0, 1]])
    q = sigma_a**2
    Q = q * np.array([[dt**4/4, 0,       dt**3/2, 0],
                      [0,       dt**4/4, 0,       dt**3/2],
                      [dt**3/2, 0,       dt**2,   0],
                      [0,       dt**3/2, 0,       dt**2]])
    H = np.array([[1, 0, 0, 0],
                  [0, 1, 0, 0]])
    R = np.diag([sigma_meas**2, sigma_meas**2])
    return F, Q, H, R

def kf_run(F, Q, H, R, x0, P0, z_seq):
    x_hat, P = x0.copy(), P0.copy()
    X = [x_hat.copy()]
    for z in z_seq:
        # Predict
        x_pred = F @ x_hat
        P_pred = F @ P @ F.T + Q
        # Update
        y = z - H @ x_pred
        S = H @ P_pred @ H.T + R
        K = P_pred @ H.T @ np.linalg.inv(S)
        x_hat = x_pred + K @ y
        P = (np.eye(len(P)) - K @ H) @ P_pred
        X.append(x_hat.copy())
    return np.array(X)

# -------------------
# Demo for SCOUT
# -------------------
if __name__ == "__main__":
    dt = 2.0
    T = 20  # 20 s 
    N = int(T/dt)

    # True trajectory: SCOUT going at 1.5 m/s east, slight north drift
    v_true = np.array([1.5, 0.2])
    truth = []
    pos = np.array([0., 0.])
    for _ in range(N+1):
        truth.append(pos.copy())
        pos += v_true * dt
    truth = np.array(truth)

    # Simulate LBL measurements with 1 m sigma
    rng = np.random.default_rng(42)
    z_meas = truth[1:] + rng.normal(0, 1.0, size=truth[1:].shape)

    # KF setup for SCOUT
    sigma_a = 0.15   # process accel noise (m/s²)
    sigma_meas = 1.0 # LBL measurement noise (m)
    F, Q, H, R = build_cv_matrices(dt, sigma_a, sigma_meas)

    # Init
    x0 = np.array([0, 0, 0, 0])  # pos, vel
    P0 = np.diag([5, 5, 1, 1])   # uncertain start

    # Run filter
    est = kf_run(F, Q, H, R, x0, P0, z_meas)

    # Compare
    plt.figure(figsize=(7,6))
    plt.plot(truth[:,0], truth[:,1], 'k-', label="True path")
    plt.scatter(z_meas[:,0], z_meas[:,1], s=20, c='g', alpha=0.5, label="LBL meas")
    plt.plot(est[:,0], est[:,1], 'b-o', label="KF estimate")
    plt.axis('equal'); plt.grid(True)
    plt.xlabel("X [m]"); plt.ylabel("Y [m]")
    plt.title("SCOUT AUV KF with noisy LBL (σ=1 m)")
    plt.legend()
    plt.show()
