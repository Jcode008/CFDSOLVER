import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# Setup paths
script_dir = Path(__file__).resolve().parent
project_dir = script_dir.parent
build_dir = project_dir / 'build'

print(f"Looking for data files in: {build_dir}")

# Get final timestep files
u_files = sorted(build_dir.glob('u_*.csv'), key=lambda p: int(p.stem.split('_')[1]))
v_files = sorted(build_dir.glob('v_*.csv'), key=lambda p: int(p.stem.split('_')[1]))
p_files = sorted(build_dir.glob('p_*.csv'), key=lambda p: int(p.stem.split('_')[1]))

if not u_files or not v_files or not p_files:
    print("Missing data files! Make sure u_*, v_*, and p_* CSV files exist.")
    sys.exit(1)

print(f"Found {len(u_files)} timesteps")

# Load the last timestep for analysis
u = np.genfromtxt(u_files[-1], delimiter=',', dtype=float)
v = np.genfromtxt(v_files[-1], delimiter=',', dtype=float)
p = np.genfromtxt(p_files[-1], delimiter=',', dtype=float)

print(f"Data shape: {u.shape}")

# Calculate velocity magnitude and dynamic pressure
vel_mag = np.sqrt(u**2 + v**2)

# Fluid properties (match solver)
rho = 1.0
U_infty = 5.0
nu = 2e-4

# Pressure coefficient: Cp = (p - p_infty) / (0.5 * rho * U_infty^2)
# Assuming p_infty ≈ mean pressure far from airfoil
p_infty = np.nanmean(p[:, -50:])  # Average pressure in outlet region
q_infty = 0.5 * rho * U_infty**2
Cp = (p - p_infty) / q_infty

print(f"\nPressure Statistics:")
print(f"  p_infty = {p_infty:.3f} Pa")
print(f"  p_min = {np.nanmin(p):.3f} Pa")
print(f"  p_max = {np.nanmax(p):.3f} Pa")
print(f"  Cp_min = {np.nanmin(Cp):.3f}")
print(f"  Cp_max = {np.nanmax(Cp):.3f}")

print(f"\nVelocity Statistics:")
print(f"  u_min = {np.nanmin(u):.3f} m/s")
print(f"  u_max = {np.nanmax(u):.3f} m/s")
print(f"  v_min = {np.nanmin(v):.3f} m/s")
print(f"  v_max = {np.nanmax(v):.3f} m/s")
print(f"  |V|_max = {np.nanmax(vel_mag):.3f} m/s")

# Create comprehensive plots
fig = plt.figure(figsize=(20, 12))

# 1. U velocity
ax1 = plt.subplot(2, 3, 1)
im1 = ax1.imshow(u, origin='lower', cmap='RdBu_r', aspect='auto', 
                 vmin=0, vmax=8, interpolation='bilinear')
plt.colorbar(im1, ax=ax1, label='u (m/s)')
ax1.set_title('U-Velocity', fontsize=14, weight='bold')
ax1.set_xlabel('Grid X')
ax1.set_ylabel('Grid Y')

# 2. V velocity
ax2 = plt.subplot(2, 3, 2)
im2 = ax2.imshow(v, origin='lower', cmap='RdBu_r', aspect='auto',
                 vmin=-2, vmax=2, interpolation='bilinear')
plt.colorbar(im2, ax=ax2, label='v (m/s)')
ax2.set_title('V-Velocity', fontsize=14, weight='bold')
ax2.set_xlabel('Grid X')
ax2.set_ylabel('Grid Y')

# 3. Velocity magnitude
ax3 = plt.subplot(2, 3, 3)
im3 = ax3.imshow(vel_mag, origin='lower', cmap='turbo', aspect='auto',
                 vmin=0, vmax=8, interpolation='bilinear')
plt.colorbar(im3, ax=ax3, label='|V| (m/s)')
ax3.set_title('Velocity Magnitude', fontsize=14, weight='bold')
ax3.set_xlabel('Grid X')
ax3.set_ylabel('Grid Y')

# 4. Pressure field
ax4 = plt.subplot(2, 3, 4)
im4 = ax4.imshow(p, origin='lower', cmap='RdBu_r', aspect='auto', interpolation='bilinear')
plt.colorbar(im4, ax=ax4, label='p (Pa)')
ax4.set_title('Pressure Field', fontsize=14, weight='bold')
ax4.set_xlabel('Grid X')
ax4.set_ylabel('Grid Y')

# 5. Pressure coefficient
ax5 = plt.subplot(2, 3, 5)
im5 = ax5.imshow(Cp, origin='lower', cmap='RdBu_r', aspect='auto',
                 vmin=-3, vmax=1, interpolation='bilinear')
plt.colorbar(im5, ax=ax5, label='Cp')
ax5.set_title('Pressure Coefficient', fontsize=14, weight='bold')
ax5.set_xlabel('Grid X')
ax5.set_ylabel('Grid Y')

# 6. Streamwise velocity comparison
ax6 = plt.subplot(2, 3, 6)
# Extract centerline velocity (y=100, middle of domain)
centerline_idx = u.shape[0] // 2
u_centerline = u[centerline_idx, :]
x_grid = np.arange(len(u_centerline))
ax6.plot(x_grid, u_centerline, 'b-', linewidth=2, label='u velocity')
ax6.axhline(U_infty, color='r', linestyle='--', linewidth=1, label=f'U∞ = {U_infty} m/s')
ax6.set_xlabel('Grid X', fontsize=12)
ax6.set_ylabel('u (m/s)', fontsize=12)
ax6.set_title('Centerline U-Velocity', fontsize=14, weight='bold')
ax6.grid(True, alpha=0.3)
ax6.legend()
ax6.set_ylim([0, 10])

plt.suptitle(f'Aerodynamic Analysis - Timestep {u_files[-1].stem.split("_")[1]}', 
             fontsize=16, weight='bold')
plt.tight_layout()

output_path = script_dir / 'aerodynamic_analysis.png'
plt.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"\nSaved: {output_path}")

# Additional analysis: Velocity profiles at different x-stations
fig2, axes = plt.subplots(2, 2, figsize=(14, 10))

# Airfoil is around x=60-120 (grid coords), analyze at different stations
stations = [50, 90, 130, 200]  # Before, mid, after, far wake
station_labels = ['Upstream', 'Mid-airfoil', 'Wake', 'Far wake']

for idx, (x_pos, label) in enumerate(zip(stations, station_labels)):
    ax = axes[idx // 2, idx % 2]
    
    if x_pos < u.shape[1]:
        u_profile = u[:, x_pos]
        v_profile = v[:, x_pos]
        y_grid = np.arange(len(u_profile))
        
        ax.plot(u_profile, y_grid, 'b-', linewidth=2, label='u velocity')
        ax.plot(v_profile, y_grid, 'r-', linewidth=2, label='v velocity')
        ax.axvline(U_infty, color='b', linestyle='--', alpha=0.5, label=f'U∞')
        ax.axvline(0, color='k', linestyle='-', alpha=0.3, linewidth=0.5)
        ax.set_xlabel('Velocity (m/s)', fontsize=12)
        ax.set_ylabel('Grid Y', fontsize=12)
        ax.set_title(f'{label} (x={x_pos})', fontsize=13, weight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend()
        ax.set_xlim([-2, 10])

plt.suptitle('Velocity Profiles at Different Stations', fontsize=15, weight='bold')
plt.tight_layout()

profile_path = script_dir / 'velocity_profiles.png'
plt.savefig(profile_path, dpi=300, bbox_inches='tight')
print(f"Saved: {profile_path}")

print("\n=== Analysis complete ===")
plt.show()
