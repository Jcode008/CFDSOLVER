import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

# Chaotic three-body problem with random initial conditions
# This creates unpredictable, beautiful motion

class ThreeBody:
    def __init__(self):
        # Masses
        self.m = np.array([1.0, 1.2, 0.9])
        
        # Initial positions (slightly chaotic configuration)
        self.pos = np.array([
            [1.0, 0.0, 0.5],
            [-0.5, 0.866, -0.3],
            [-0.5, -0.866, 0.2]
        ])
        
        # Initial velocities (give them a kick for chaos!)
        self.vel = np.array([
            [0.2, 0.5, 0.3],
            [-0.4, 0.1, -0.2],
            [0.3, -0.6, 0.4]
        ])
        
        self.G = 1.0  # Gravitational constant
        self.dt = 0.01  # Time step
        
    def forces(self):
        """Calculate gravitational forces between all bodies"""
        F = np.zeros((3, 3))
        
        for i in range(3):
            for j in range(3):
                if i != j:
                    r_vec = self.pos[j] - self.pos[i]
                    r = np.linalg.norm(r_vec)
                    # Softening to prevent singularities
                    r = max(r, 0.1)
                    F[i] += self.G * self.m[i] * self.m[j] * r_vec / (r**3)
        
        return F
    
    def step(self):
        """Leapfrog integration for stability"""
        # Calculate forces
        F = self.forces()
        
        # Update velocities (half step)
        acc = F / self.m[:, np.newaxis]
        self.vel += 0.5 * acc * self.dt
        
        # Update positions
        self.pos += self.vel * self.dt
        
        # Recalculate forces
        F = self.forces()
        acc = F / self.m[:, np.newaxis]
        
        # Update velocities (second half step)
        self.vel += 0.5 * acc * self.dt

# Simulation parameters
fps = 30
duration = 20.0  # 20 seconds of chaos!
n_frames = int(fps * duration)
trail_length = 150  # Longer trails for chaotic beauty

# Run simulation
print("Simulating chaotic orbits...")
sim = ThreeBody()
positions = np.zeros((n_frames, 3, 3))

for i in range(n_frames):
    positions[i] = sim.pos.copy()
    sim.step()
    if i % 100 == 0:
        print(f"  Frame {i}/{n_frames}")

print("Creating 3D animation...")

# Setup the figure with custom background
fig = plt.figure(figsize=(12, 12), facecolor='#282a36')
ax = fig.add_subplot(111, projection='3d', facecolor='#282a36')

# Set viewing angle
ax.view_init(elev=20, azim=45)

# Find bounds for nice scaling
all_pos = positions.reshape(-1, 3)
margin = 0.5
x_range = [all_pos[:, 0].min() - margin, all_pos[:, 0].max() + margin]
y_range = [all_pos[:, 1].min() - margin, all_pos[:, 1].max() + margin]
z_range = [all_pos[:, 2].min() - margin, all_pos[:, 2].max() + margin]

ax.set_xlim(x_range)
ax.set_ylim(y_range)
ax.set_zlim(z_range)

# Dark theme styling
ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False
ax.xaxis.pane.set_edgecolor('#44475a')
ax.yaxis.pane.set_edgecolor('#44475a')
ax.zaxis.pane.set_edgecolor('#44475a')
ax.grid(True, alpha=0.2, color='#6272a4')

# Hide axes labels for cleaner look
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_zticklabels([])

# Colors for the three bodies (Dracula theme inspired)
colors = ['#ff79c6', '#8be9fd', '#50fa7b']  # Pink, Cyan, Green
planet_names = ['Crimson', 'Azure', 'Emerald']

# Create planet scatter plots
planets = [ax.scatter([], [], [], c=color, s=400, alpha=0.9, 
                     edgecolors='#f8f8f2', linewidth=2, marker='o')
          for color in colors]

# Create glowing halos
halos = [ax.scatter([], [], [], c=color, s=800, alpha=0.15, marker='o')
         for color in colors]

# Create trails (fading paths)
trails = [ax.plot([], [], [], color=color, linewidth=2, alpha=0.7)[0] 
         for color in colors]

# Faint full trails showing complete history
full_trails = [ax.plot([], [], [], color=color, linewidth=0.5, alpha=0.2)[0] 
              for color in colors]

# Time counter
time_text = fig.text(0.5, 0.05, '', ha='center', va='center', 
                    fontsize=12, color='#bd93f9', family='monospace')

def init():
    for planet in planets:
        planet._offsets3d = ([], [], [])
    for halo in halos:
        halo._offsets3d = ([], [], [])
    for trail in trails:
        trail.set_data([], [])
        trail.set_3d_properties([])
    return planets + halos + trails

def animate(frame):
    # Rotate camera for dynamic view
    ax.view_init(elev=20 + 10*np.sin(frame/50), 
                azim=45 + frame * 0.5)
    
    # Update planets
    for i, (planet, halo) in enumerate(zip(planets, halos)):
        x, y, z = positions[frame, i]
        planet._offsets3d = ([x], [y], [z])
        halo._offsets3d = ([x], [y], [z])
    
    # Update recent trails (fading)
    start = max(0, frame - trail_length)
    for i, trail in enumerate(trails):
        trail_x = positions[start:frame+1, i, 0]
        trail_y = positions[start:frame+1, i, 1]
        trail_z = positions[start:frame+1, i, 2]
        trail.set_data(trail_x, trail_y)
        trail.set_3d_properties(trail_z)
    
    # Update full trails (entire history, very faint)
    for i, full_trail in enumerate(full_trails):
        trail_x = positions[:frame+1, i, 0]
        trail_y = positions[:frame+1, i, 1]
        trail_z = positions[:frame+1, i, 2]
        full_trail.set_data(trail_x, trail_y)
        full_trail.set_3d_properties(trail_z)
    
    # Update time counter
    time_text.set_text(f'T = {frame * sim.dt:.2f}s')
    
    # Pulse halos
    pulse = 0.1 + 0.1 * np.sin(2 * np.pi * frame / 30)
    for halo in halos:
        halo.set_alpha(pulse)
    
    return planets + halos + trails + full_trails + [time_text]

# Create animation
print("Rendering frames...")
anim = animation.FuncAnimation(fig, animate, init_func=init,
                              frames=n_frames, interval=1000/fps, 
                              blit=False, repeat=True)

# Save as GIF
print("Saving GIF (this will take a few minutes due to 3D complexity)...")
writer = animation.PillowWriter(fps=fps, metadata=dict(artist='CFD Solver'), 
                               bitrate=2000)
anim.save('three_body_3d.gif', writer=writer, dpi=80)

print("✓ Saved as three_body_3d.gif")
print(f"  Duration: {duration}s")
print(f"  Frames: {n_frames}")
print(f"  FPS: {fps}")
print(f"  Background: #282a36 (Dracula theme)")

plt.close()
print("\n🌌 Chaotic 3D animation complete!")
