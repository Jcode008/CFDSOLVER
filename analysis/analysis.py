import sys
import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# Establish project locations relative to this script
script_dir = Path(__file__).resolve().parent
project_dir = script_dir.parent
candidate_dirs = [
    project_dir / 'build',  # Check build directory FIRST
    # project_dir / 'build' / 'Debug',
    project_dir / 'snapshots',
    project_dir
]

build_dir = next((d for d in candidate_dirs if d.exists()), None)
if build_dir is None:
    print("No build directory with snapshot files found. Re-run the solver first.")
    sys.exit(1)

print(f"Looking for CSV files in: {build_dir}")

# Get all u_*.csv files sorted by timestep number
files = sorted(build_dir.glob('u_*.csv'), key=lambda path: int(path.stem.split('_')[1]))
print(f"Found {len(files)} timestep files")

if not files:
    print("No snapshot files found. Run the solver first or check the snapshot path.")
    sys.exit(0)

# Load all timesteps
timesteps = []
timestep_numbers = []
reference_shape = None
for file in files:
    data = np.genfromtxt(
        file,
        delimiter=',',
        dtype=float,
        filling_values=np.nan,
        invalid_raise=False
    )
    if data.ndim == 1:
        data = data[np.newaxis, :]

    if reference_shape is None:
        reference_shape = data.shape
    elif data.shape != reference_shape:
        print(f"Skipping {file} (shape {data.shape} vs expected {reference_shape})")
        continue

    timesteps.append(data)
    timestep_num = int(file.stem.split('_')[1])
    timestep_numbers.append(timestep_num)

# Convert to 3D array: [time, y, x]
u_all = np.array(timesteps)
print(f"Data shape: {u_all.shape}")

# Check for NaNs
print("\n=== NaN Check ===")
for i, data in enumerate(timesteps):
    if np.isnan(data).any():
        print(f"⚠️  NaNs found at timestep {timestep_numbers[i]}: {np.isnan(data).sum()} cells")
    else:
        print(f"✓ Timestep {timestep_numbers[i]}: No NaNs")

# Calculate statistics over time
print("\n=== Statistics ===")
mins = [np.nanmin(data) for data in timesteps]
maxs = [np.nanmax(data) for data in timesteps]
means = [np.nanmean(data) for data in timesteps]

# Plot 1: Min/Max/Mean over time
plt.figure(figsize=(12, 5))
plt.plot(timestep_numbers, mins, label='min', linewidth=2)
plt.plot(timestep_numbers, maxs, label='max', linewidth=2)
plt.plot(timestep_numbers, means, label='mean', linewidth=2, linestyle='--')
plt.xlabel('Timestep', fontsize=12)
plt.ylabel('u velocity (m/s)', fontsize=12)
plt.title('U-Velocity Statistics Over Time', fontsize=14)
plt.legend(fontsize=10)
plt.grid(True, alpha=0.3)
plt.tight_layout()
output_dir = script_dir
output_dir.mkdir(parents=True, exist_ok=True)
stats_path = output_dir / 'velocity_statistics.png'
plt.savefig(stats_path, dpi=300)
print(f"Saved: {stats_path}")

# Define visualization range - focus on velocity deficit in wake
# Center around freestream velocity to see acceleration AND deceleration
vmin, vmax = 2.0, 10.0  # Wider range to reduce saturation flickering, focus on wake
cmap = 'turbo'  # Better for seeing flow structures

# Plot 2: First timestep (initial condition)
plt.figure(figsize=(16, 7))
im = plt.imshow(u_all[0], origin='lower', cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax, interpolation='bilinear')
cbar = plt.colorbar(im, label='u velocity (m/s)', extend='both')
cbar.ax.tick_params(labelsize=10)
plt.title(f'U-Velocity Field - Timestep {timestep_numbers[0]}', fontsize=16, weight='bold')
plt.xlabel('Grid X', fontsize=14)
plt.ylabel('Grid Y', fontsize=14)
plt.tight_layout()
initial_path = output_dir / 'velocity_field_initial.png'
plt.savefig(initial_path, dpi=300, bbox_inches='tight')
print(f"Saved: {initial_path}")

# Plot 3: Last timestep (final state) with enhanced contrast
plt.figure(figsize=(16, 7))
im = plt.imshow(u_all[-1], origin='lower', cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax, interpolation='bilinear')
cbar = plt.colorbar(im, label='u velocity (m/s)', extend='both')
cbar.ax.tick_params(labelsize=10)
plt.title(f'U-Velocity Field - Timestep {timestep_numbers[-1]}', fontsize=16, weight='bold')
plt.xlabel('Grid X', fontsize=14)
plt.ylabel('Grid Y', fontsize=14)
plt.tight_layout()
final_path = output_dir / 'velocity_field_final.png'
plt.savefig(final_path, dpi=300, bbox_inches='tight')
print(f"Saved: {final_path}")

# Plot 4: Contour plot with more levels for detail
plt.figure(figsize=(16, 7))
levels = np.linspace(vmin, vmax, 50)  # More levels = smoother gradients
cf = plt.contourf(u_all[-1], levels=levels, cmap=cmap, extend='both')
cbar = plt.colorbar(cf, label='u velocity (m/s)')
cbar.ax.tick_params(labelsize=10)
# Add contour lines for structure
plt.contour(u_all[-1], levels=15, colors='black', alpha=0.3, linewidths=0.5)
plt.title(f'U-Velocity Contours - Timestep {timestep_numbers[-1]}', fontsize=16, weight='bold')
plt.xlabel('Grid X', fontsize=14)
plt.ylabel('Grid Y', fontsize=14)
plt.tight_layout()
contour_path = output_dir / 'velocity_contours.png'
plt.savefig(contour_path, dpi=300, bbox_inches='tight')
print(f"Saved: {contour_path}")

# Plot 5: Multiple timesteps comparison with consistent scaling
fig, axes = plt.subplots(2, 3, figsize=(20, 11))
axes = axes.flatten()

# Select 6 evenly spaced timesteps
indices = np.linspace(0, len(timesteps)-1, 6, dtype=int)

for idx, ax in enumerate(axes):
    i = indices[idx]
    im = ax.imshow(u_all[i], origin='lower', cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax, interpolation='bilinear')
    ax.set_title(f'Timestep {timestep_numbers[i]}', fontsize=14, weight='bold')
    ax.set_xlabel('Grid X', fontsize=12)
    ax.set_ylabel('Grid Y', fontsize=12)
    cbar = plt.colorbar(im, ax=ax, label='u (m/s)')
    cbar.ax.tick_params(labelsize=9)

plt.suptitle('U-Velocity Evolution', fontsize=18, weight='bold')
plt.tight_layout()
evolution_path = output_dir / 'velocity_evolution.png'
plt.savefig(evolution_path, dpi=300, bbox_inches='tight')
print(f"Saved: {evolution_path}")

print("\n=== All plots saved to analysis/ directory ===")

# ===========================
# ANIMATION - Create MP4 video
# ===========================
print("\n=== Creating Animation ===")
import matplotlib.animation as animation

# Set up the figure and axis for animation
fig_anim, ax_anim = plt.subplots(figsize=(16, 8))

# Create initial plot
im_anim = ax_anim.imshow(u_all[0], origin='lower', cmap=cmap, aspect='auto', 
                         vmin=vmin, vmax=vmax, interpolation='bilinear')
cbar_anim = plt.colorbar(im_anim, ax=ax_anim, label='u velocity (m/s)', extend='both')
cbar_anim.ax.tick_params(labelsize=10)

title_anim = ax_anim.set_title(f'Flow Evolution - Timestep {timestep_numbers[0]}', 
                               fontsize=16, weight='bold')
ax_anim.set_xlabel('Grid X', fontsize=14)
ax_anim.set_ylabel('Grid Y', fontsize=14)

# Animation update function
def update(frame):
    im_anim.set_data(u_all[frame])
    title_anim.set_text(f'Flow Evolution - Timestep {timestep_numbers[frame]}')
    return [im_anim, title_anim]

# Create animation
print(f"Generating animation with {len(timesteps)} frames...")
anim = animation.FuncAnimation(
    fig_anim, 
    update, 
    frames=len(timesteps),
    interval=50,  # 50ms between frames = 20 fps
    blit=True
)

# Save as MP4
animation_path = output_dir / 'flow_animation.mp4'
print(f"Saving animation to {animation_path}...")
print("This might take a minute...")

try:
    # Try ffmpeg first
    try:
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=20, metadata=dict(artist='CFD Solver'), bitrate=3000)
        anim.save(str(animation_path), writer=writer, dpi=150)
        print(f"✅ Animation saved: {animation_path}")
    except:
        # Fallback to pillow (saves as GIF instead)
        print("ffmpeg not available, saving as GIF instead...")
        gif_path = output_dir / 'flow_animation.gif'
        anim.save(str(gif_path), writer='pillow', fps=20, dpi=100)
        print(f"✅ Animation saved as GIF: {gif_path}")
        animation_path = gif_path
    
    # Try to open the video/gif
    import os
    try:
        os.startfile(str(animation_path))
        print("Opening animation...")
    except:
        print(f"Could not auto-open. Manually open: {animation_path}")
        
except Exception as e:
    print(f"⚠️  Could not create animation: {e}")
    import traceback
    traceback.print_exc()

plt.close(fig_anim)

print("\n=== Opening static plots ===")
plt.show()


