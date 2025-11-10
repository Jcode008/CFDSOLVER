"""
Create a beautiful diagram showing O-grid topology advantage
Highlights: no wake cut, radial structure, smooth wrapping
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

# Load O-grid from build directory
build_dir = os.path.join(os.path.dirname(__file__), '..', 'build')
grid_path = os.path.join(build_dir, 'grid.csv')

print("Loading O-grid...")
grid = pd.read_csv(grid_path)

nxi = grid['i'].max() + 1
neta = grid['j'].max() + 1
print(f"Grid size: {nxi} x {neta}")

# Reshape into 2D arrays
x = grid['x'].values.reshape(neta, nxi)
y = grid['y'].values.reshape(neta, nxi)

# Create figure
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# ============================================================
# LEFT: Full O-grid with radial structure highlighted
# ============================================================
ax1 = axes[0]
ax1.set_aspect('equal')
ax1.set_title('O-Grid Topology: Wraps 0→2π (No Wake Cut)', fontsize=16, fontweight='bold')

# Plot radial circles (constant η)
for j in range(0, neta, max(1, neta//10)):
    alpha_val = 0.3 + 0.4 * (j / neta)  # Fade out toward farfield
    linewidth = 1.0 if j < 10 else 0.5
    ax1.plot(x[j, :], y[j, :], 'b-', alpha=alpha_val, linewidth=linewidth)

# Plot radial spokes (constant ξ) - highlight key angles
key_angles = [0, nxi//4, nxi//2, 3*nxi//4]  # 0°, 90°, 180°, 270°
for i in key_angles:
    ax1.plot(x[:, i], y[:, i], 'r-', alpha=0.8, linewidth=1.5)

# Plot some intermediate spokes
for i in range(0, nxi, max(1, nxi//24)):
    if i not in key_angles:
        ax1.plot(x[:, i], y[:, i], 'r-', alpha=0.3, linewidth=0.5)

# Highlight airfoil
ax1.plot(x[0, :], y[0, :], 'k-', linewidth=3, label='NACA 2412', zorder=10)

# Add annotations
ax1.annotate('Leading Edge', xy=(x[0, nxi//2], y[0, nxi//2]), 
             xytext=(x[0, nxi//2]-0.3, y[0, nxi//2]-0.5),
             fontsize=12, fontweight='bold',
             arrowprops=dict(arrowstyle='->', lw=2, color='green'))

ax1.annotate('Trailing Edge', xy=(x[0, 0], y[0, 0]), 
             xytext=(x[0, 0]+0.3, y[0, 0]+0.3),
             fontsize=12, fontweight='bold',
             arrowprops=dict(arrowstyle='->', lw=2, color='green'))

# Add text box explaining topology
textstr = 'ξ: circumferential (wraps 0→2π)\nη: radial (surface→farfield)\n\n✓ No discontinuity at wake\n✓ Smooth radial lines\n✓ Better aspect ratio'
props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
ax1.text(0.02, 0.98, textstr, transform=ax1.transAxes, fontsize=11,
         verticalalignment='top', bbox=props, family='monospace')

ax1.grid(True, alpha=0.3)
ax1.legend(loc='lower right', fontsize=12)
ax1.set_xlabel('x (m)', fontsize=12)
ax1.set_ylabel('y (m)', fontsize=12)

# ============================================================
# RIGHT: Zoom near airfoil showing boundary layer clustering
# ============================================================
ax2 = axes[1]
ax2.set_aspect('equal')
ax2.set_title('Near-Airfoil Clustering (First 15 Radial Layers)', fontsize=16, fontweight='bold')

# Plot first 15 radial layers with different colors
n_layers = min(15, neta)
colors = plt.cm.viridis(np.linspace(0, 1, n_layers))

for j in range(n_layers):
    ax2.plot(x[j, :], y[j, :], color=colors[j], linewidth=1.5, alpha=0.8)

# Plot radial spokes
for i in range(0, nxi, max(1, nxi//40)):
    ax2.plot(x[:n_layers, i], y[:n_layers, i], 'r-', alpha=0.4, linewidth=0.5)

# Highlight airfoil surface
ax2.plot(x[0, :], y[0, :], 'k-', linewidth=3, label='Airfoil surface', zorder=10)

# Show first cell height
i_mid = nxi // 2
dx = x[1, i_mid] - x[0, i_mid]
dy = y[1, i_mid] - y[0, i_mid]
first_cell = np.sqrt(dx**2 + dy**2)

ax2.plot([x[0, i_mid], x[1, i_mid]], [y[0, i_mid], y[1, i_mid]], 
         'g-', linewidth=3, label=f'1st cell: {first_cell*1000:.2f} mm')

# Add info box
textstr = f'Radial clustering (β=4.0)\n\nFirst cell: {first_cell*1000:.2f} mm\nSmooth growth to farfield\n\nIdeal for boundary layer!'
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.8)
ax2.text(0.98, 0.02, textstr, transform=ax2.transAxes, fontsize=11,
         verticalalignment='bottom', horizontalalignment='right', 
         bbox=props, family='monospace')

ax2.grid(True, alpha=0.3)
ax2.legend(loc='upper left', fontsize=11)
ax2.set_xlabel('x (m)', fontsize=12)
ax2.set_ylabel('y (m)', fontsize=12)
ax2.set_xlim(-0.15, 0.75)
ax2.set_ylim(-0.35, 0.35)

# Save figure
plt.tight_layout()
output_path = os.path.join(build_dir, 'ogrid_topology.png')
plt.savefig(output_path, dpi=200, bbox_inches='tight')
print(f"\n✓ O-grid topology diagram saved to {output_path}")

print("\nO-Grid Statistics:")
print("=" * 50)
print(f"  Resolution: {nxi} × {neta} = {nxi*neta:,} points")
print(f"  Farfield: {np.max(np.sqrt(x**2 + y**2)):.2f} m")
print(f"  First cell height: {first_cell*1000:.2f} mm")
print(f"  Topology: Wraps 0→2π around airfoil")
print(f"  Wake cut: NONE! (smooth periodic boundary)")
