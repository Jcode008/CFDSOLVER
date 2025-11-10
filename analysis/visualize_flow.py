"""
Visualize CFD solution on curvilinear O-grid
Shows velocity field, pressure, streamlines around NACA 2412
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

# Load grid and solution from build directory
build_dir = os.path.join(os.path.dirname(__file__), '..', 'build')

print("Loading grid and solution...")
grid_df = pd.read_csv(os.path.join(build_dir, 'grid.csv'))
u_df = pd.read_csv(os.path.join(build_dir, 'u_final.csv'))
v_df = pd.read_csv(os.path.join(build_dir, 'v_final.csv'))
p_df = pd.read_csv(os.path.join(build_dir, 'p_final.csv'))

# Get grid dimensions
nxi = grid_df['i'].max() + 1
neta = grid_df['j'].max() + 1
print(f"Grid: {nxi} × {neta} = {nxi*neta:,} points")

# Reshape arrays
x = grid_df['x'].values.reshape(neta, nxi)
y = grid_df['y'].values.reshape(neta, nxi)
u = u_df['value'].values.reshape(neta, nxi)
v = v_df['value'].values.reshape(neta, nxi)
p = p_df['value'].values.reshape(neta, nxi)

# Compute derived quantities
velocity_mag = np.sqrt(u**2 + v**2)
vorticity = np.zeros_like(u)

# Simple vorticity computation (ω_z = ∂v/∂x - ∂u/∂y)
for j in range(1, neta-1):
    for i in range(1, nxi-1):
        dvdx = (v[j, i+1] - v[j, i-1]) / (x[j, i+1] - x[j, i-1])
        dudy = (u[j+1, i] - u[j-1, i]) / (y[j+1, i] - y[j-1, i])
        vorticity[j, i] = dvdx - dudy

print(f"Velocity range: [{velocity_mag.min():.2f}, {velocity_mag.max():.2f}] m/s")
print(f"Pressure range: [{p.min():.2f}, {p.max():.2f}] Pa")

# Create comprehensive figure
fig = plt.figure(figsize=(18, 12))

# ============================================================
# Plot 1: Velocity magnitude with streamlines
# ============================================================
ax1 = plt.subplot(2, 3, 1)
levels = np.linspace(0, velocity_mag.max(), 30)
cf1 = ax1.contourf(x, y, velocity_mag, levels=levels, cmap='jet', extend='max')
ax1.contour(x, y, velocity_mag, levels=10, colors='k', alpha=0.3, linewidths=0.5)

# Add velocity vectors instead of streamlines (works on curvilinear grids)
skip_i = max(1, nxi // 20)
skip_j = max(1, neta // 10)
i_indices = np.arange(0, nxi, skip_i)
j_indices = np.arange(0, neta, skip_j)
X_sub = x[np.ix_(j_indices, i_indices)]
Y_sub = y[np.ix_(j_indices, i_indices)]
U_sub = u[np.ix_(j_indices, i_indices)]
V_sub = v[np.ix_(j_indices, i_indices)]
ax1.quiver(X_sub, Y_sub, U_sub, V_sub, color='white', 
           scale=60, width=0.002, alpha=0.7)

ax1.plot(x[0, :], y[0, :], 'k-', linewidth=2.5)  # Airfoil
ax1.set_aspect('equal')
ax1.set_title('Velocity Magnitude + Velocity Vectors', fontsize=14, fontweight='bold')
ax1.set_xlabel('x (m)')
ax1.set_ylabel('y (m)')
ax1.set_xlim(-0.3, 1.0)
ax1.set_ylim(-0.6, 0.6)
plt.colorbar(cf1, ax=ax1, label='|V| (m/s)')

# ============================================================
# Plot 2: Pressure coefficient
# ============================================================
ax2 = plt.subplot(2, 3, 2)
U_inf = 5.0  # Freestream velocity
rho = 1.0    # Density (assumed)
q_inf = 0.5 * rho * U_inf**2
Cp = p / q_inf

levels_cp = np.linspace(-2, 1.5, 30)
cf2 = ax2.contourf(x, y, Cp, levels=levels_cp, cmap='RdBu_r', extend='both')
ax2.contour(x, y, Cp, levels=15, colors='k', alpha=0.3, linewidths=0.5)
ax2.plot(x[0, :], y[0, :], 'k-', linewidth=2.5)
ax2.set_aspect('equal')
ax2.set_title('Pressure Coefficient (Cp)', fontsize=14, fontweight='bold')
ax2.set_xlabel('x (m)')
ax2.set_ylabel('y (m)')
ax2.set_xlim(-0.3, 1.0)
ax2.set_ylim(-0.6, 0.6)
plt.colorbar(cf2, ax=ax2, label='Cp = p/q∞')

# ============================================================
# Plot 3: Vorticity (shows boundary layer and wake)
# ============================================================
ax3 = plt.subplot(2, 3, 3)
vort_max = np.percentile(np.abs(vorticity), 98)  # Clip extreme values
levels_vort = np.linspace(-vort_max, vort_max, 30)
cf3 = ax3.contourf(x, y, vorticity, levels=levels_vort, cmap='RdBu_r', extend='both')
ax3.plot(x[0, :], y[0, :], 'k-', linewidth=2.5)
ax3.set_aspect('equal')
ax3.set_title('Vorticity (ωz = ∂v/∂x - ∂u/∂y)', fontsize=14, fontweight='bold')
ax3.set_xlabel('x (m)')
ax3.set_ylabel('y (m)')
ax3.set_xlim(-0.3, 1.0)
ax3.set_ylim(-0.6, 0.6)
plt.colorbar(cf3, ax=ax3, label='ωz (1/s)')

# ============================================================
# Plot 4: Cp distribution on airfoil surface
# ============================================================
ax4 = plt.subplot(2, 3, 4)
x_surface = x[0, :]
Cp_surface = Cp[0, :]

# Split into upper and lower surfaces
mid_idx = nxi // 2
x_lower = x_surface[:mid_idx+1]
x_upper = x_surface[mid_idx:]
Cp_lower = Cp_surface[:mid_idx+1]
Cp_upper = Cp_surface[mid_idx:]

ax4.plot(x_lower / 0.6, Cp_lower, 'b-o', linewidth=2, markersize=3, label='Lower surface')
ax4.plot(x_upper / 0.6, Cp_upper, 'r-s', linewidth=2, markersize=3, label='Upper surface')
ax4.axhline(0, color='k', linestyle='--', alpha=0.3)
ax4.set_xlabel('x/c', fontsize=12)
ax4.set_ylabel('Cp', fontsize=12)
ax4.set_title('Surface Pressure Distribution', fontsize=14, fontweight='bold')
ax4.legend()
ax4.grid(True, alpha=0.3)
ax4.invert_yaxis()  # Convention: negative Cp upward

# ============================================================
# Plot 5: Velocity profiles at different chordwise stations
# ============================================================
ax5 = plt.subplot(2, 3, 5)

# Plot velocity profiles at x/c = 0.25, 0.5, 0.75
stations = [0.25, 0.5, 0.75]
colors = ['blue', 'green', 'red']

for station, color in zip(stations, colors):
    # Find closest i index
    target_x = station * 0.6
    i_station = np.argmin(np.abs(x[0, :] - target_x))
    
    # Extract velocity profile
    y_profile = y[:, i_station]
    u_profile = u[:, i_station]
    
    # Normalize: distance from surface, velocity by freestream
    y_norm = y_profile - y[0, i_station]
    u_norm = u_profile / U_inf
    
    # Plot first 20 points (boundary layer region)
    ax5.plot(u_norm[:20], y_norm[:20]*1000, '-o', color=color, 
             linewidth=2, markersize=4, label=f'x/c = {station:.2f}')

ax5.set_xlabel('u / U∞', fontsize=12)
ax5.set_ylabel('Distance from surface (mm)', fontsize=12)
ax5.set_title('Boundary Layer Velocity Profiles', fontsize=14, fontweight='bold')
ax5.legend()
ax5.grid(True, alpha=0.3)
ax5.set_xlim(0, 1.2)

# ============================================================
# Plot 6: Velocity vectors (quiver plot near airfoil)
# ============================================================
ax6 = plt.subplot(2, 3, 6)

# Subsample and plot velocity vectors
skip_i = max(1, nxi // 30)
skip_j = 2
i_indices = np.arange(0, nxi, skip_i)
j_indices = np.arange(0, min(25, neta), skip_j)

X_sub = x[np.ix_(j_indices, i_indices)]
Y_sub = y[np.ix_(j_indices, i_indices)]
U_sub = u[np.ix_(j_indices, i_indices)]
V_sub = v[np.ix_(j_indices, i_indices)]

# Color by velocity magnitude
vel_mag_sub = np.sqrt(U_sub**2 + V_sub**2)
quiv = ax6.quiver(X_sub, Y_sub, U_sub, V_sub, vel_mag_sub,
                   cmap='jet', scale=60, width=0.003, alpha=0.8)
ax6.plot(x[0, :], y[0, :], 'k-', linewidth=2.5)
ax6.set_aspect('equal')
ax6.set_title('Velocity Vectors (Near Airfoil)', fontsize=14, fontweight='bold')
ax6.set_xlabel('x (m)')
ax6.set_ylabel('y (m)')
ax6.set_xlim(-0.1, 0.7)
ax6.set_ylim(-0.3, 0.3)
plt.colorbar(quiv, ax=ax6, label='|V| (m/s)')

# ============================================================
# Save and display
# ============================================================
plt.tight_layout()
output_path = os.path.join(build_dir, 'flow_solution.png')
plt.savefig(output_path, dpi=200, bbox_inches='tight')
print(f"\n✓ Flow visualization saved to {output_path}")

# Print summary statistics
print("\n" + "="*60)
print("FLOW SOLUTION SUMMARY (NACA 2412 at α=-5°, Re=15000)")
print("="*60)
print(f"Velocity:  min={velocity_mag.min():.3f} m/s, max={velocity_mag.max():.3f} m/s")
print(f"Pressure:  min={p.min():.2f} Pa, max={p.max():.2f} Pa")
print(f"Cp:        min={Cp.min():.3f}, max={Cp.max():.3f}")
print(f"Vorticity: min={vorticity.min():.2f} 1/s, max={vorticity.max():.2f} 1/s")

# Estimate forces (simple integration)
# Integrate pressure on surface
chord = 0.6
dx_surface = np.diff(x[0, :])
dy_surface = np.diff(y[0, :])
ds_surface = np.sqrt(dx_surface**2 + dy_surface**2)
p_avg = (Cp[0, :-1] + Cp[0, 1:]) / 2.0

# Normal vectors (pointing into fluid)
nx = -dy_surface / ds_surface
ny = dx_surface / ds_surface

# Force components (pressure only, simplified)
Fx = np.sum(p_avg * nx * ds_surface) * q_inf
Fy = np.sum(p_avg * ny * ds_surface) * q_inf

# Rotate to lift/drag (α = -5°)
alpha = -5.0 * np.pi / 180.0
L = -Fx * np.sin(alpha) + Fy * np.cos(alpha)
D = Fx * np.cos(alpha) + Fy * np.sin(alpha)

# Coefficients
S = chord * 1.0  # Unit span
CL = L / (q_inf * S)
CD = D / (q_inf * S)

print(f"\nEstimated Aerodynamic Coefficients:")
print(f"  CL = {CL:.4f}")
print(f"  CD = {CD:.4f}")
print(f"  L/D = {CL/CD if CD != 0 else 0:.2f}")
print("\n(Note: Pressure forces only - add viscous drag for full CD)")
print("="*60)
