"""
Analyze results from Cartesian grid simulation
Creates comprehensive visualization and computes aerodynamic coefficients
"""
import numpy as np
import matplotlib.pyplot as plt
import os

# Setup
build_dir = r'C:\Users\graha\CFDSolver\build\Release'
print("Analyzing Cartesian grid simulation results...\n")

# Simulation parameters
Lx, Ly = 4.0, 2.0
nx, ny = 800, 400
rho = 1.225  # kg/m³
U_infty = 51.44  # m/s (100 knots)
chord = 0.6  # m
alpha = -5.0  # degrees (in code) = +5° pitch up

# Create grid
x = np.linspace(0, Lx, nx)
y = np.linspace(0, Ly, ny)
X, Y = np.meshgrid(x, y)
dx = Lx / (nx - 1)
dy = Ly / (ny - 1)

# Load final timestep
final_step = 4950
print(f"Loading timestep {final_step}...")
u = np.genfromtxt(os.path.join(build_dir, f'u_{final_step}.csv'), delimiter=',')
v = np.genfromtxt(os.path.join(build_dir, f'v_{final_step}.csv'), delimiter=',')
p = np.genfromtxt(os.path.join(build_dir, f'p_{final_step}.csv'), delimiter=',')

print(f"Grid: {nx} × {ny} = {nx*ny:,} cells")
print(f"Data shape: u={u.shape}, v={v.shape}, p={p.shape}\n")

# Compute derived quantities
vel_mag = np.sqrt(u**2 + v**2)
q_infty = 0.5 * rho * U_infty**2  # Dynamic pressure
Cp = p / q_infty  # Pressure coefficient

# Mask for valid (fluid) cells
valid = ~np.isnan(u)

# Print statistics
print("=" * 60)
print("FLOW STATISTICS")
print("=" * 60)
print(f"Velocity:  min={np.nanmin(vel_mag):.2f} m/s, max={np.nanmax(vel_mag):.2f} m/s")
print(f"Pressure:  min={np.nanmin(p):.2f} Pa, max={np.nanmax(p):.2f} Pa")
print(f"Cp:        min={np.nanmin(Cp):.3f}, max={np.nanmax(Cp):.3f}")
print(f"Dynamic pressure: {q_infty:.2f} Pa")
print()

# Create comprehensive visualization
fig = plt.figure(figsize=(18, 12))

# 1. Velocity magnitude
ax1 = plt.subplot(2, 3, 1)
levels_v = np.linspace(0, np.nanmax(vel_mag), 30)
cf1 = ax1.contourf(X, Y, vel_mag, levels=levels_v, cmap='jet', extend='max')
ax1.set_xlabel('x (m)')
ax1.set_ylabel('y (m)')
ax1.set_title('Velocity Magnitude |V| (m/s)', fontweight='bold')
ax1.set_aspect('equal')
plt.colorbar(cf1, ax=ax1, label='|V| (m/s)')
ax1.set_xlim(0.2, 1.8)
ax1.set_ylim(0.5, 1.5)

# 2. Pressure field
ax2 = plt.subplot(2, 3, 2)
levels_p = np.linspace(np.nanmin(p), np.nanmax(p), 30)
cf2 = ax2.contourf(X, Y, p, levels=levels_p, cmap='RdBu_r', extend='both')
ax2.set_xlabel('x (m)')
ax2.set_ylabel('y (m)')
ax2.set_title('Pressure (Pa)', fontweight='bold')
ax2.set_aspect('equal')
plt.colorbar(cf2, ax=ax2, label='p (Pa)')
ax2.set_xlim(0.2, 1.8)
ax2.set_ylim(0.5, 1.5)

# 3. Pressure coefficient
ax3 = plt.subplot(2, 3, 3)
levels_cp = np.linspace(-2, 1, 30)
cf3 = ax3.contourf(X, Y, Cp, levels=levels_cp, cmap='RdBu_r', extend='both')
ax3.set_xlabel('x (m)')
ax3.set_ylabel('y (m)')
ax3.set_title('Pressure Coefficient Cp', fontweight='bold')
ax3.set_aspect('equal')
plt.colorbar(cf3, ax=ax3, label='Cp')
ax3.set_xlim(0.2, 1.8)
ax3.set_ylim(0.5, 1.5)

# 4. Streamlines
ax4 = plt.subplot(2, 3, 4)
# Subsample for clearer streamlines
step = 4
ax4.streamplot(X[::step, ::step], Y[::step, ::step], 
               u[::step, ::step], v[::step, ::step],
               color=vel_mag[::step, ::step], cmap='viridis',
               linewidth=1, density=1.5, arrowsize=1)
ax4.set_xlabel('x (m)')
ax4.set_ylabel('y (m)')
ax4.set_title('Streamlines', fontweight='bold')
ax4.set_aspect('equal')
ax4.set_xlim(0.2, 1.8)
ax4.set_ylim(0.5, 1.5)

# 5. Velocity vectors (zoomed near airfoil)
ax5 = plt.subplot(2, 3, 5)
step = 10
X_zoom = X[100:300, 300:600]
Y_zoom = Y[100:300, 300:600]
u_zoom = u[100:300, 300:600]
v_zoom = v[100:300, 300:600]
vel_zoom = vel_mag[100:300, 300:600]
ax5.quiver(X_zoom[::step, ::step], Y_zoom[::step, ::step],
           u_zoom[::step, ::step], v_zoom[::step, ::step],
           vel_zoom[::step, ::step], cmap='jet', scale=500)
ax5.set_xlabel('x (m)')
ax5.set_ylabel('y (m)')
ax5.set_title('Velocity Vectors (Near Airfoil)', fontweight='bold')
ax5.set_aspect('equal')

# 6. Vorticity
ax6 = plt.subplot(2, 3, 6)
# Compute vorticity: ω = dv/dx - du/dy
vorticity = np.zeros_like(u)
for i in range(1, ny-1):
    for j in range(1, nx-1):
        if valid[i, j]:
            dvdx = (v[i, j+1] - v[i, j-1]) / (2*dx)
            dudy = (u[i+1, j] - u[i-1, j]) / (2*dy)
            vorticity[i, j] = dvdx - dudy

vorticity[~valid] = np.nan
vort_lim = np.nanpercentile(np.abs(vorticity), 99)  # Use 99th percentile
levels_w = np.linspace(-vort_lim, vort_lim, 30)
cf6 = ax6.contourf(X, Y, vorticity, levels=levels_w, cmap='RdBu_r', extend='both')
ax6.set_xlabel('x (m)')
ax6.set_ylabel('y (m)')
ax6.set_title('Vorticity ω (1/s)', fontweight='bold')
ax6.set_aspect('equal')
plt.colorbar(cf6, ax=ax6, label='ω (1/s)')
ax6.set_xlim(0.2, 1.8)
ax6.set_ylim(0.5, 1.5)

# Add simulation info
info_text = f'''NACA 2412 Airfoil @ 100 knots
α = +5° (nose up), Re = 2.08×10⁶
Grid: {nx}×{ny}, Time: {final_step*1.25e-5:.4f}s
ISA Sea Level Conditions'''

fig.text(0.5, 0.96, info_text, ha='center', fontsize=11,
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

plt.tight_layout(rect=[0, 0, 1, 0.95])

# Save figure
output_path = os.path.join(build_dir, 'complete_analysis_100kts.png')
plt.savefig(output_path, dpi=200, bbox_inches='tight')
print(f"✓ Complete analysis saved to: {output_path}\n")

# Simple force estimation (very rough)
print("=" * 60)
print("ROUGH AERODYNAMIC ESTIMATE")
print("=" * 60)
print("Note: This is a 2D simulation with simplified boundary conditions")
print("Results are qualitative, not quantitative!")
print()

# Find airfoil region (where pressure is defined but velocity is zero)
airfoil_mask = np.isnan(u)

# Try to estimate forces (very approximate)
if np.any(airfoil_mask):
    # Count airfoil cells
    n_airfoil = np.sum(airfoil_mask)
    print(f"Airfoil cells detected: {n_airfoil}")
    
    # Average pressure on airfoil boundary
    # This is super crude - just for qualitative understanding
    boundary_pressures = []
    for i in range(1, ny-1):
        for j in range(1, nx-1):
            if airfoil_mask[i, j]:
                # Check if it's a boundary cell (has fluid neighbor)
                if (not airfoil_mask[i-1, j] or not airfoil_mask[i+1, j] or 
                    not airfoil_mask[i, j-1] or not airfoil_mask[i, j+1]):
                    # Average pressure from fluid neighbors
                    p_neighbors = []
                    if not airfoil_mask[i-1, j]: p_neighbors.append(p[i-1, j])
                    if not airfoil_mask[i+1, j]: p_neighbors.append(p[i+1, j])
                    if not airfoil_mask[i, j-1]: p_neighbors.append(p[i, j-1])
                    if not airfoil_mask[i, j+1]: p_neighbors.append(p[i, j+1])
                    if p_neighbors:
                        boundary_pressures.append(np.mean(p_neighbors))
    
    if boundary_pressures:
        p_avg = np.mean(boundary_pressures)
        print(f"Average boundary pressure: {p_avg:.2f} Pa")
        print(f"Average Cp on airfoil: {p_avg/q_infty:.3f}")

print("\n" + "=" * 60)
print("✓ Analysis complete!")
print("=" * 60)
