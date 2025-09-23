#!/usr/bin/env python3
import math

# --- Inputs ---
# Final AUV position from Part A.1 (meters)
auv_pos = (8.84633979220194, 3.74017247479654)

# LBL transponder positions (meters)
transponders = {
    "T1 (100,-100)": (100, -100),
    "T2 (-100,100)": (-100, 100),
    "T3 (100,100)":  (100, 100),
    "T4 (-100,-100)":(-100, -100),
}

c = 1500.0  # speed of sound in water [m/s]
round_trip = False  # set True if you want two-way TOF

# --- Compute distances and TOF ---
print(f"AUV final position: x={auv_pos[0]:.4f} m, y={auv_pos[1]:.4f} m\n")
print(f"{'Transponder':18s}  {'Range [m]':>10s}  {'TOF [s]':>10s}  {'TOF [ms]':>10s}")
for name, (bx, by) in transponders.items():
    dx = auv_pos[0] - bx
    dy = auv_pos[1] - by
    r = math.hypot(dx, dy)                # geometric range [m]
    tof = (2.0 * r / c) if round_trip else (r / c)  # time of flight
    print(f"{name:18s}  {r:10.3f}  {tof:10.6f}  {1e3*tof:10.3f}")
