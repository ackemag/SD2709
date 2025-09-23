#!/usr/bin/python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def wrap_pi(a):
    return np.arctan2(np.sin(a), np.cos(a))

# ---- Unicycle step with (v, omega) as inputs ----
def unicycle_motion(x, y, theta, v, omega, dt):
    """
    Inputs
      x,y,theta : current pose (m, m, rad) in world frame
      v, omega  : surge speed (m/s) and yaw rate (rad/s), given in body frame
      dt        : step [s]
    Outputs
      x_new, y_new, theta_new : pose after dt
    """
    x_new = x + v * np.cos(theta) * dt
    y_new = y + v * np.sin(theta) * dt
    theta_new = wrap_pi(theta + omega * dt)
    return x_new, y_new, theta_new



# ---- Initial conditions & controls ----
x, y, theta = 0.0, 0.0, 0.0       # start at origin, heading +x
v, omega = 1.0, 0.1               # m/s, rad/s (can be arrays per step)
dt, T = 2.0, 10.0                 # s
steps = int(T / dt)

trajectory = [(0.0, 0.0, 0.0)]
for k in range(steps):
    x, y, theta = unicycle_motion(x, y, theta, v, omega, dt)
    trajectory.append((x, y, theta))


# ---- Print a quick table ----
df = pd.DataFrame(trajectory, columns=["x [m]", "y [m]", "theta [rad]"])
df.index = (df.index * dt)  # time stamps
print(df)

# ---- Plot ----
plt.figure(figsize=(7,5))
plt.plot(df["x [m]"], df["y [m]"], marker='o', label='DR path')
plt.scatter(df["x [m]"].iloc[0], df["y [m]"].iloc[0], label='Start', zorder=5)
plt.scatter(df["x [m]"].iloc[-1], df["y [m]"].iloc[-1], label='End', zorder=5)
plt.title('AUV Dead-Reckoning (2D Unicycle)')
plt.xlabel('X [m]'); plt.ylabel('Y [m]')
plt.grid(True); plt.legend()
plt.show()
