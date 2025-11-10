"""
Create animation of flow evolution over time
Shows velocity magnitude developing around airfoil
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd
import os

# Setup
build_dir = os.path.join(os.path.dirname(__file__), '..', 'build')
print("Creating flow evolution animation...")

# Load grid
grid = pd.read_csv(os.path.join(build_dir, 'grid.csv'))
nxi = grid['i'].max() + 1
neta = grid['j'].max() + 1
x = grid['x'].values.reshape(neta, nxi)
y = grid['y'].values.reshape(neta, nxi)

# Find timesteps
timesteps = []
for f in os.listdir(build_dir):
    if f.startswith('u_') and f.endswith('.csv') and f != 'u_final.csv':
        try:
            step = int(f[2:-4])
            timesteps.append(step)
        except:
            pass
timesteps.sort()

if len(timesteps) == 0:
    print("No timestep data found! Run solver with intermediate exports.")
    exit(1)

print(f"Found {len(timesteps)} timesteps: {timesteps[0]} to {timesteps[-1]}")

# Create figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

# Setup plot 1: Velocity magnitude
ax1.set_aspect('equal')
ax1.set_xlabel('x (m)', fontsize=12)
ax1.set_ylabel('y (m)', fontsize=12)
ax1.set_xlim(-0.2, 1.0)
ax1.set_ylim(-0.5, 0.5)
ax1.plot(x[0, :], y[0, :], 'k-', linewidth=3, zorder=10)

# Setup plot 2: Pressure field
ax2.set_aspect('equal')
ax2.set_xlabel('x (m)', fontsize=12)
ax2.set_ylabel('y (m)', fontsize=12)
ax2.set_xlim(-0.2, 1.0)
ax2.set_ylim(-0.5, 0.5)
ax2.plot(x[0, :], y[0, :], 'k-', linewidth=3, zorder=10)

# Initialize with first timestep
step = timesteps[0]
u = pd.read_csv(os.path.join(build_dir, f'u_{step}.csv'))['value'].values.reshape(neta, nxi)
v = pd.read_csv(os.path.join(build_dir, f'v_{step}.csv'))['value'].values.reshape(neta, nxi)
p = pd.read_csv(os.path.join(build_dir, f'p_{step}.csv'))['value'].values.reshape(neta, nxi)
vel_mag = np.sqrt(u**2 + v**2)

# Find global min/max for consistent color scales
print("Computing global ranges...")
v_max = 0
p_min, p_max = 0, 0
for step in timesteps[::5]:  # Sample every 5th
    u_t = pd.read_csv(os.path.join(build_dir, f'u_{step}.csv'))['value'].values.reshape(neta, nxi)
    v_t = pd.read_csv(os.path.join(build_dir, f'v_{step}.csv'))['value'].values.reshape(neta, nxi)
    p_t = pd.read_csv(os.path.join(build_dir, f'p_{step}.csv'))['value'].values.reshape(neta, nxi)
    v_mag_t = np.sqrt(u_t**2 + v_t**2)
    v_max = max(v_max, v_mag_t.max())
    p_min = min(p_min, p_t.min())
    p_max = max(p_max, p_t.max())

print(f"Velocity range: [0, {v_max:.2f}] m/s")
print(f"Pressure range: [{p_min:.2f}, {p_max:.2f}] Pa")

levels_v = np.linspace(0, v_max, 25)
levels_p = np.linspace(p_min, p_max, 25)

cf1 = ax1.contourf(x, y, vel_mag, levels=levels_v, cmap='jet', extend='max')
cf2 = ax2.contourf(x, y, p, levels=levels_p, cmap='RdBu_r', extend='both')

cbar1 = plt.colorbar(cf1, ax=ax1, label='|V| (m/s)')
cbar2 = plt.colorbar(cf2, ax=ax2, label='p (Pa)')

title1 = ax1.set_title(f'Velocity Magnitude - Step {step}', fontsize=14, fontweight='bold')
title2 = ax2.set_title(f'Pressure Field - Step {step}', fontsize=14, fontweight='bold')

def animate(frame):
    """Update animation frame"""
    step = timesteps[frame]
    
    # Load data
    u = pd.read_csv(os.path.join(build_dir, f'u_{step}.csv'))['value'].values.reshape(neta, nxi)
    v = pd.read_csv(os.path.join(build_dir, f'v_{step}.csv'))['value'].values.reshape(neta, nxi)
    p = pd.read_csv(os.path.join(build_dir, f'p_{step}.csv'))['value'].values.reshape(neta, nxi)
    vel_mag = np.sqrt(u**2 + v**2)
    
    # Clear and redraw
    for ax in [ax1, ax2]:
        for coll in ax.collections:
            coll.remove()
    
    ax1.contourf(x, y, vel_mag, levels=levels_v, cmap='jet', extend='max')
    ax2.contourf(x, y, p, levels=levels_p, cmap='RdBu_r', extend='both')
    
    # Update titles
    title1.set_text(f'Velocity Magnitude - Step {step} (t={step*1e-5:.4f}s)')
    title2.set_text(f'Pressure Field - Step {step} (t={step*1e-5:.4f}s)')
    
    # Re-draw airfoil
    ax1.plot(x[0, :], y[0, :], 'k-', linewidth=3, zorder=10)
    ax2.plot(x[0, :], y[0, :], 'k-', linewidth=3, zorder=10)
    
    print(f"  Frame {frame+1}/{len(timesteps)}: step {step}")
    
    return [title1, title2]

# Create animation
print(f"\nGenerating animation with {len(timesteps)} frames...")
anim = animation.FuncAnimation(fig, animate, frames=len(timesteps), 
                               interval=200, blit=True, repeat=True)

# Save as GIF
output_path = os.path.join(build_dir, 'flow_evolution.gif')
print(f"\nSaving animation to {output_path}...")
anim.save(output_path, writer='pillow', fps=5, dpi=100)

file_size = os.path.getsize(output_path) / 1024 / 1024
print(f"\n✓ Animation saved! ({file_size:.1f} MB)")
print(f"  Frames: {len(timesteps)}")
print(f"  Duration: {len(timesteps)/5:.1f} seconds @ 5 fps")
