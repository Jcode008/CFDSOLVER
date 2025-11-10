import sys
import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Establish project locations relative to this script
script_dir = Path(__file__).resolve().parent
project_dir = script_dir.parent
build_dir = project_dir / 'build'

print(f"Looking for CSV files in: {build_dir}")

# Get all u_*.csv files sorted by timestep number
files = sorted(build_dir.glob('u_*.csv'), key=lambda path: int(path.stem.split('_')[1]))
print(f"Found {len(files)} timestep files")

if not files:
    print("No snapshot files found. Run the solver first!")
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
print(f"Creating animation with {len(timesteps)} frames...")

# Set up the figure and axis
fig, ax = plt.subplots(figsize=(16, 8))

# Define visualization range
vmin, vmax = 0.0, 8.0
cmap = 'turbo'

# Create initial plot
im = ax.imshow(u_all[0], origin='lower', cmap=cmap, aspect='auto', 
               vmin=vmin, vmax=vmax, interpolation='bilinear')
cbar = plt.colorbar(im, ax=ax, label='u velocity (m/s)', extend='both')
cbar.ax.tick_params(labelsize=10)

title = ax.set_title(f'Flow Evolution - Timestep {timestep_numbers[0]}', 
                     fontsize=16, weight='bold')
ax.set_xlabel('Grid X', fontsize=14)
ax.set_ylabel('Grid Y', fontsize=14)

# Animation update function
def update(frame):
    im.set_data(u_all[frame])
    title.set_text(f'Flow Evolution - Timestep {timestep_numbers[frame]}')
    return [im, title]

# Create animation
print("Generating animation...")
anim = animation.FuncAnimation(
    fig, 
    update, 
    frames=len(timesteps),
    interval=50,  # 50ms between frames = 20 fps
    blit=True
)

# Save as MP4
output_path = script_dir / 'flow_animation.mp4'
print(f"Saving animation to {output_path}...")
print("This might take a minute...")

# Use ffmpeg writer
Writer = animation.writers['ffmpeg']
writer = Writer(fps=20, metadata=dict(artist='CFD Solver'), bitrate=3000)

anim.save(str(output_path), writer=writer, dpi=150)

print(f"\n✅ Animation saved to: {output_path}")
print("Opening video...")

# Try to open the video
try:
    os.startfile(str(output_path))
except:
    print(f"Could not auto-open. Manually open: {output_path}")

plt.close()
