"""
Visualize body-fitted O-grid
Shows grid lines and quality metrics for O-grid topology
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

# Load grid from build directory
build_dir = os.path.join(os.path.dirname(__file__), '..', 'build')
grid_path = os.path.join(build_dir, 'grid.csv')
metrics_path = os.path.join(build_dir, 'metrics.csv')

print("Loading grid...")
grid = pd.read_csv(grid_path)
metrics = pd.read_csv(metrics_path)

nxi = grid['i'].max() + 1
neta = grid['j'].max() + 1
print(f"Grid size: {nxi} x {neta}")

# Reshape into 2D arrays
x = grid['x'].values.reshape(neta, nxi)
y = grid['y'].values.reshape(neta, nxi)

# Create figure with subplots
fig = plt.figure(figsize=(16, 10))

# ============================================================
# Plot 1: Full O-grid
# ============================================================
ax1 = plt.subplot(2, 3, 1)
ax1.set_aspect('equal')
ax1.set_title('Full O-Grid (No Wake Cut!)', fontsize=14, fontweight='bold')

# Plot grid lines in ξ direction (around airfoil)
for j in range(0, neta, max(1, neta//15)):  # Every radial layer
    ax1.plot(x[j, :], y[j, :], 'b-', alpha=0.5, linewidth=0.5)

# Plot grid lines in η direction (radial lines)
for i in range(0, nxi, max(1, nxi//30)):  # Radial spokes
    ax1.plot(x[:, i], y[:, i], 'r-', alpha=0.5, linewidth=0.5)

# Highlight airfoil (η=0)
ax1.plot(x[0, :], y[0, :], 'k-', linewidth=2.5, label='NACA 2412 airfoil')
ax1.grid(True, alpha=0.3)
ax1.legend()
ax1.set_xlabel('x (m)')
ax1.set_ylabel('y (m)')

# ============================================================
# Plot 2: Zoom on airfoil (radial structure)
# ============================================================
ax2 = plt.subplot(2, 3, 2)
ax2.set_aspect('equal')
ax2.set_title('O-Grid Near Airfoil (Radial Lines)', fontsize=14, fontweight='bold')

# Plot every grid line near airfoil
for j in range(0, min(20, neta), 2):  # First 20 radial layers
    ax2.plot(x[j, :], y[j, :], 'b-', alpha=0.6, linewidth=0.5)

for i in range(0, nxi, max(1, nxi//60)):  # Radial spokes
    ax2.plot(x[:20, i], y[:20, i], 'r-', alpha=0.6, linewidth=0.5)

ax2.plot(x[0, :], y[0, :], 'k-', linewidth=2.5, label='NACA 2412')
ax2.grid(True, alpha=0.3)
ax2.legend()
ax2.set_xlabel('x (m)')
ax2.set_ylabel('y (m)')
ax2.set_xlim(-0.1, 0.7)
ax2.set_ylim(-0.3, 0.3)

# ============================================================
# Plot 3: Jacobian field
# ============================================================
ax3 = plt.subplot(2, 3, 3)
J = metrics['J'].values.reshape(neta, nxi)
im = ax3.contourf(x, y, J, levels=50, cmap='viridis')
ax3.plot(x[0, :], y[0, :], 'k-', linewidth=2)
ax3.set_aspect('equal')
ax3.set_title('Jacobian Field', fontsize=14, fontweight='bold')
ax3.set_xlabel('x (m)')
ax3.set_ylabel('y (m)')
plt.colorbar(im, ax=ax3, label='J = x_ξ·y_η - x_η·y_ξ')
ax3.set_xlim(-1, 2)
ax3.set_ylim(-1.5, 1.5)

# Check for negative Jacobians
if np.any(J <= 0):
    print(f"⚠ WARNING: {np.sum(J <= 0)} points have J ≤ 0 (grid folded!)")
else:
    print(f"✓ All Jacobians positive (min={J.min():.6f}, max={J.max():.6f})")

# ============================================================
# Plot 4: Grid orthogonality
# ============================================================
ax4 = plt.subplot(2, 3, 4)

# Compute angle between grid lines
# ∇ξ · ∇η = 0 for orthogonal grid
xi_x = metrics['xi_x'].values.reshape(neta, nxi)
xi_y = metrics['xi_y'].values.reshape(neta, nxi)
eta_x = metrics['eta_x'].values.reshape(neta, nxi)
eta_y = metrics['eta_y'].values.reshape(neta, nxi)

# Dot product of gradients
dot_product = xi_x * eta_x + xi_y * eta_y
# Magnitudes
mag_xi = np.sqrt(xi_x**2 + xi_y**2)
mag_eta = np.sqrt(eta_x**2 + eta_y**2)
# Angle (cos θ = dot / (|a||b|))
cos_theta = dot_product / (mag_xi * mag_eta + 1e-10)
angle_deg = np.abs(90.0 - np.degrees(np.arccos(np.clip(cos_theta, -1, 1))))

im = ax4.contourf(x, y, angle_deg, levels=np.linspace(0, 30, 31), cmap='RdYlGn_r')
ax4.plot(x[0, :], y[0, :], 'k-', linewidth=2)
ax4.set_aspect('equal')
ax4.set_title('Grid Orthogonality Error', fontsize=14, fontweight='bold')
ax4.set_xlabel('x (m)')
ax4.set_ylabel('y (m)')
plt.colorbar(im, ax=ax4, label='Deviation from 90° (degrees)')
ax4.set_xlim(-0.2, 1.0)
ax4.set_ylim(-0.5, 0.5)

print(f"Orthogonality: mean error = {np.mean(angle_deg):.2f}°, max = {np.max(angle_deg):.2f}°")

# ============================================================
# Plot 5: Grid spacing (η direction)
# ============================================================
ax5 = plt.subplot(2, 3, 5)

# Compute spacing in η direction at midpoint of airfoil
i_mid = nxi // 2
spacing = []
eta_coords = []
for j in range(neta - 1):
    dx = x[j+1, i_mid] - x[j, i_mid]
    dy = y[j+1, i_mid] - y[j, i_mid]
    ds = np.sqrt(dx**2 + dy**2)
    spacing.append(ds)
    eta_coords.append(j / (neta - 1))

ax5.semilogy(eta_coords, spacing, 'b-o', linewidth=2, markersize=4)
ax5.grid(True, alpha=0.3)
ax5.set_title('Normal Spacing (at mid-chord)', fontsize=14, fontweight='bold')
ax5.set_xlabel('η (0=surface, 1=farfield)')
ax5.set_ylabel('Grid spacing Δs (m)')

# Check first cell height (important for boundary layer)
first_cell = spacing[0]
print(f"First cell height (wall): {first_cell*1000:.3f} mm")

# ============================================================
# Plot 6: Grid aspect ratio
# ============================================================
ax6 = plt.subplot(2, 3, 6)

# Compute cell aspect ratio (Δξ / Δη)
aspect_ratios = []
for j in range(1, min(30, neta-1)):  # Near wall
    for i in range(1, nxi-1):
        # ξ direction length
        dx_xi = x[j, i+1] - x[j, i-1]
        dy_xi = y[j, i+1] - y[j, i-1]
        ds_xi = np.sqrt(dx_xi**2 + dy_xi**2) / 2.0
        
        # η direction length
        dx_eta = x[j+1, i] - x[j-1, i]
        dy_eta = y[j+1, i] - y[j-1, i]
        ds_eta = np.sqrt(dx_eta**2 + dy_eta**2) / 2.0
        
        ar = ds_xi / (ds_eta + 1e-10)
        if ar < 100:  # Filter outliers
            aspect_ratios.append(ar)

ax6.hist(aspect_ratios, bins=50, edgecolor='black', alpha=0.7)
ax6.axvline(1.0, color='r', linestyle='--', linewidth=2, label='Ideal (AR=1)')
ax6.set_title('Cell Aspect Ratio Distribution', fontsize=14, fontweight='bold')
ax6.set_xlabel('Aspect Ratio (Δξ / Δη)')
ax6.set_ylabel('Number of cells')
ax6.legend()
ax6.grid(True, alpha=0.3)

median_ar = np.median(aspect_ratios)
print(f"Cell aspect ratio: median = {median_ar:.2f}")

# ============================================================
# Summary
# ============================================================
plt.tight_layout()
output_path = os.path.join(build_dir, 'grid_analysis.png')
plt.savefig(output_path, dpi=150, bbox_inches='tight')
print(f"\n✓ Grid visualization saved to {output_path}")
# plt.show()  # Skip interactive display

print("\nGrid Quality Summary:")
print("=" * 50)
print(f"  Total points: {nxi * neta:,}")
print(f"  Jacobian: [{J.min():.6f}, {J.max():.6f}]")
print(f"  Orthogonality error: {np.mean(angle_deg):.2f}° ± {np.std(angle_deg):.2f}°")
print(f"  First cell height: {first_cell*1000:.3f} mm")
print(f"  Median aspect ratio: {median_ar:.2f}")

if J.min() > 0 and np.mean(angle_deg) < 15 and median_ar < 50:
    print("\n✓ Grid quality looks good! Ready for CFD.")
else:
    print("\n⚠ Grid may need adjustment:")
    if J.min() <= 0:
        print("  - Negative Jacobian detected")
    if np.mean(angle_deg) > 15:
        print("  - Poor orthogonality (add control functions)")
    if median_ar > 50:
        print("  - High aspect ratio (adjust clustering)")
