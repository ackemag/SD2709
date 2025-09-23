#!/usr/bin/env python3
# Kalman Filter (EKF) for a 2D unicycle with LBL (x,y) measurements
# Chatgpt was used in helping write this code
import numpy as np
import matplotlib.pyplot as plt

def wrap_pi(a):
    return np.arctan2(np.sin(a), np.cos(a))

# ---------- Motion model ----------
def f_unicycle(x, u, dt):
    """ Nonlinear state propagation for unicycle: x=[x,y,theta], u=[v,omega] """
    xk, yk, th = x
    v, om = u
    xkp1 = xk + v*np.cos(th)*dt
    ykp1 = yk + v*np.sin(th)*dt
    thp1 = wrap_pi(th + om*dt)
    return np.array([xkp1, ykp1, thp1])

def F_jacobian(x, u, dt):
    """ df/dx at (x,u) """
    _, _, th = x
    v, _ = u
    return np.array([
        [1.0, 0.0, -v*np.sin(th)*dt],
        [0.0, 1.0,  v*np.cos(th)*dt],
        [0.0, 0.0,  1.0]
    ])

def B_jacobian(x, dt):
    """ df/du at (x,u) """
    _, _, th = x
    return np.array([
        [np.cos(th)*dt, 0.0],
        [np.sin(th)*dt, 0.0],
        [0.0,           dt ]
    ])

# ---------- EKF steps ----------
def ekf_predict(x, P, u, Q, dt):
    F = F_jacobian(x, u, dt)
    B = B_jacobian(x, dt)
    x_pred = f_unicycle(x, u, dt)
    # Inject process noise through controls (speed & yaw-rate uncertainty)
    P_pred = F @ P @ F.T + B @ Q @ B.T
    return x_pred, P_pred

def ekf_update(x_pred, P_pred, z, H, R):
    # Innovation
    y = z - (H @ x_pred)
    S = H @ P_pred @ H.T + R
    K = P_pred @ H.T @ np.linalg.inv(S)
    # State/covariance update (Joseph form for numerical stability)
    x_upd = x_pred + K @ y
    I = np.eye(P_pred.shape[0])
    P_upd = (I - K @ H) @ P_pred @ (I - K @ H).T + K @ R @ K.T
    # Keep angle wrapped
    x_upd[2] = wrap_pi(x_upd[2])
    return x_upd, P_upd

# ---------- Demo / minimal sim ----------
if __name__ == "__main__":
    # Timing & inputs (same spirit as Part A.1)
    dt = 2.0
    T  = 10.0
    N  = int(T/dt)
    u  = np.array([1.0, 0.1])  # v=1 m/s, omega=0.1 rad/s (can make time-varying if you want)

    # True trajectory (for generating synthetic LBL fixes)
    x_true = np.array([0.0, 0.0, 0.0])
    truth = [x_true.copy()]
    for _ in range(N):
        x_true = f_unicycle(x_true, u, dt)
        truth.append(x_true.copy())
    truth = np.array(truth)

    # Simulate LBL: z = [x,y] + noise   (one fix per step here)
    H = np.array([[1.0, 0.0, 0.0],
                  [0.0, 1.0, 0.0]])
    # Measurement std (tune): e.g. 0.8 m 1σ in x,y
    R = np.diag([0.8**2, 0.8**2])

    rng = np.random.default_rng(7)
    lbl_meas = []
    for k in range(1, N+1):  # measurements aligned with steps 1..N
        z = truth[k, :2] + rng.normal(0.0, [np.sqrt(R[0,0]), np.sqrt(R[1,1])])
        lbl_meas.append(z)
    lbl_meas = np.array(lbl_meas)

    # EKF initialization
    x_hat = np.array([0.0, 0.0, 0.0])  # start at origin
    P     = np.diag([1.0**2, 1.0**2, (10.0*np.pi/180.0)**2])  # 1 m pos, 10 deg heading 1σ

    # Process noise on inputs (speed & yaw-rate) → mapped with B
    # Example: σ_v = 0.15 m/s, σ_ω = 0.03 rad/s
    Q_u = np.diag([0.15**2, 0.03**2])

    est_kf = [x_hat.copy()]
    est_dr = [x_hat.copy()]  # pure dead-reckoning from estimate (no measurement updates)

    for k in range(N):
        # Predict
        x_pred, P_pred = ekf_predict(x_hat, P, u, Q_u, dt)
        # DR track (just for plotting/comparison)
        est_dr.append(x_pred.copy())

        # Update with LBL (x,y)
        z_k = lbl_meas[k]
        x_hat, P = ekf_update(x_pred, P_pred, z_k, H, R)

        est_kf.append(x_hat.copy())

    est_kf = np.array(est_kf)
    est_dr = np.array(est_dr)

    # ---------- Plots ----------
    plt.figure(figsize=(7,6))
    plt.plot(truth[:,0], truth[:,1], 'k-',  label='True (DR with exact inputs)')
    plt.plot(est_dr[:,0], est_dr[:,1], 'r-o', label='Dead-reckoning (predict only)', alpha=0.6)
    plt.plot(est_kf[:,0], est_kf[:,1], 'b-o', label='EKF estimate')
    plt.scatter(lbl_meas[:,0], lbl_meas[:,1], c='g', s=25, label='LBL (x,y) measurements', zorder=5)
    plt.axis('equal'); plt.grid(True)
    plt.xlabel('X [m]'); plt.ylabel('Y [m]')
    plt.title('2D Unicycle EKF with LBL (x,y) Measurements')
    plt.legend()
    plt.tight_layout()
    plt.show()

    # Print final states
    print("\nFinal states:")
    print(f"True:  x={truth[-1,0]:.3f}, y={truth[-1,1]:.3f}, th={truth[-1,2]:.3f} rad")
    print(f"DR  :  x={est_dr[-1,0]:.3f}, y={est_dr[-1,1]:.3f}, th={est_dr[-1,2]:.3f} rad")
    print(f"EKF :  x={est_kf[-1,0]:.3f}, y={est_kf[-1,1]:.3f}, th={est_kf[-1,2]:.3f} rad")
