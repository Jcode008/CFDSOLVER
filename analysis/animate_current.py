"""
Create animation from current Cartesian grid simulation
Works with u_N.csv, v_N.csv, p_N.csv format (raw CSV, no headers)
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
import glob

# Setup
build_dir = r'C:\Users\graha\CFDSolver\build\Release'
print("Creating flow animation from Cartesian grid data...")

# Find timesteps
u_files = glob.glob(os.path.join(build_dir, 'u_*.csv'))
timesteps = []
for f in u_files:
    fname = os.path.basename(f)
    if fname == 'u_final.csv':
        continue
    try:
        step = int(fname[2:-4])
        timesteps.append(step)
    except:
        pass

timesteps.sort()
print(f"Found {len(timesteps)} timesteps: {timesteps[0]} to {timesteps[-1]}")

if len(timesteps) == 0:
    print("No data files found!")
    exit(1)

# Load first frame to get dimensions
u_0 = np.genfromtxt(os.path.join(build_dir, f'u_{timesteps[0]}.csv'), delimiter=',')
ny, nx = u_0.shape
print(f"Grid size: {nx} × {ny}")

# Create coordinate arrays
Lx, Ly = 4.0, 2.0
x = np.linspace(0, Lx, nx)
y = np.linspace(0, Ly, ny)
X, Y = np.meshgrid(x, y)

# Find global ranges for consistent colors
print("Computing value ranges...")
u_max, v_max, p_min, p_max = 0, 0, 0, 0
for step in timesteps[::max(1, len(timesteps)//10)]:  # Sample 10 frames
    u = np.genfromtxt(os.path.join(build_dir, f'u_{step}.csv'), delimiter=',')
    v = np.genfromtxt(os.path.join(build_dir, f'v_{step}.csv'), delimiter=',')
    p = np.genfromtxt(os.path.join(build_dir, f'p_{step}.csv'), delimiter=',')
    
    u_valid = u[~np.isnan(u)]
    v_valid = v[~np.isnan(v)]
    p_valid = p[~np.isnan(p)]
    
    if len(u_valid) > 0:
        u_max = max(u_max, np.abs(u_valid).max())
    if len(v_valid) > 0:
        v_max = max(v_max, np.abs(v_valid).max())
    if len(p_valid) > 0:
        p_min = min(p_min, p_valid.min())
        p_max = max(p_max, p_valid.max())

vel_max = max(u_max, v_max) * 1.1
print(f"Velocity range: [0, {vel_max:.2f}] m/s")
print(f"Pressure range: [{p_min:.2f}, {p_max:.2f}] Pa")

# Create figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

# Setup axes
for ax in [ax1, ax2]:
    ax.set_aspect('equal')
    ax.set_xlabel('x (m)', fontsize=12)
    ax.set_ylabel('y (m)', fontsize=12)
    ax.set_xlim(0, Lx)
    ax.set_ylim(0, Ly)

# Load first frame
u = np.genfromtxt(os.path.join(build_dir, f'u_{timesteps[0]}.csv'), delimiter=',')
v = np.genfromtxt(os.path.join(build_dir, f'v_{timesteps[0]}.csv'), delimiter=',')
p = np.genfromtxt(os.path.join(build_dir, f'p_{timesteps[0]}.csv'), delimiter=',')
vel_mag = np.sqrt(u**2 + v**2)

# Create initial plots
levels_v = np.linspace(0, vel_max, 30)
levels_p = np.linspace(p_min, p_max, 30)

cf1 = ax1.contourf(X, Y, vel_mag, levels=levels_v, cmap='jet', extend='max')
cf2 = ax2.contourf(X, Y, p, levels=levels_p, cmap='RdBu_r', extend='both')

cbar1 = plt.colorbar(cf1, ax=ax1, label='|V| (m/s)')
cbar2 = plt.colorbar(cf2, ax=ax2, label='Pressure (Pa)')

# Titles
dt = 1.25e-5  # From your simulation
title1 = ax1.set_title(f'Velocity Magnitude - Step {timesteps[0]} (t=0.000s)', 
                       fontsize=14, fontweight='bold')
title2 = ax2.set_title(f'Pressure Field - Step {timesteps[0]} (t=0.000s)', 
                       fontsize=14, fontweight='bold')

# Add simulation info
info_text = f'NACA 2412, α=5°, U∞=51.4 m/s (100 kts)\nRe=2.1×10⁶'
fig.text(0.5, 0.02, info_text, ha='center', fontsize=10, 
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

def animate(frame):
    """Update animation frame"""
    step = timesteps[frame]
    t = step * dt
    
    # Load data
    u = np.genfromtxt(os.path.join(build_dir, f'u_{step}.csv'), delimiter=',')
    v = np.genfromtxt(os.path.join(build_dir, f'v_{step}.csv'), delimiter=',')
    p = np.genfromtxt(os.path.join(build_dir, f'p_{step}.csv'), delimiter=',')
    vel_mag = np.sqrt(u**2 + v**2)
    
    # Clear old contours
    for coll in ax1.collections:
        coll.remove()
    for coll in ax2.collections:
        coll.remove()
    
    # Redraw
    ax1.contourf(X, Y, vel_mag, levels=levels_v, cmap='jet', extend='max')
    ax2.contourf(X, Y, p, levels=levels_p, cmap='RdBu_r', extend='both')
    
    # Update titles
    title1.set_text(f'Velocity Magnitude - Step {step} (t={t:.4f}s)')
    title2.set_text(f'Pressure Field - Step {step} (t={t:.4f}s)')
    
    print(f"  Frame {frame+1}/{len(timesteps)}: step {step}")
    
    return [title1, title2]

# Create animation
print(f"\nGenerating animation with {len(timesteps)} frames...")
anim = animation.FuncAnimation(fig, animate, frames=len(timesteps), 
                               interval=150, blit=True, repeat=True)

# Save
output_path = os.path.join(build_dir, 'flow_animation_100kts.gif')
print(f"\nSaving animation to {output_path}...")
print("This may take a few minutes...")

anim.save(output_path, writer='pillow', fps=6, dpi=100)

file_size = os.path.getsize(output_path) / 1024 / 1024
print(f"\n✓ Animation complete! ({file_size:.1f} MB)")
print(f"  Frames: {len(timesteps)}")
print(f"  Duration: {len(timesteps)/6:.1f} seconds @ 6 fps")
print(f"  Location: {output_path}")
