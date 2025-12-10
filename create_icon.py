"""
Create a simple icon for the CFD Solver application
Generates a .ico file from a simple matplotlib plot
"""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrow
import io
from PIL import Image

# Create figure
fig, ax = plt.subplots(1, 1, figsize=(2, 2), dpi=128)
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis('off')

# Draw simple airfoil shape
theta = np.linspace(0, np.pi, 100)
x_upper = 0.5 + 0.35 * np.cos(theta)
y_upper = 0.5 + 0.15 * np.sin(theta) + 0.03

x_lower = 0.5 + 0.35 * np.cos(theta[::-1])
y_lower = 0.5 - 0.12 * np.sin(theta[::-1])

x_airfoil = np.concatenate([x_upper, x_lower])
y_airfoil = np.concatenate([y_upper, y_lower])

ax.fill(x_airfoil, y_airfoil, color='#2E86AB', edgecolor='#06283D', linewidth=2)

# Draw flow lines
for y_pos in [0.25, 0.5, 0.75]:
    ax.arrow(0.05, y_pos, 0.3, 0, head_width=0.03, head_length=0.05, 
             fc='#A23E48', ec='#A23E48', linewidth=1.5, alpha=0.7)

# Add title
ax.text(0.5, 0.05, 'CFD', ha='center', va='center', 
        fontsize=20, fontweight='bold', color='#06283D')

plt.tight_layout(pad=0)

# Save as PNG first
png_path = 'icon_temp.png'
plt.savefig(png_path, dpi=256, transparent=True, bbox_inches='tight', pad_inches=0)
plt.close()

# Convert to ICO
try:
    img = Image.open(png_path)
    
    # Resize to common icon sizes
    sizes = [(256, 256), (128, 128), (64, 64), (48, 48), (32, 32), (16, 16)]
    icon_images = []
    
    for size in sizes:
        icon_img = img.resize(size, Image.Resampling.LANCZOS)
        icon_images.append(icon_img)
    
    # Save as ICO
    icon_images[0].save('app_icon.ico', format='ICO', sizes=[(s[0], s[1]) for s in sizes])
    
    print("✓ Icon created: app_icon.ico")
    print("  Sizes: 16x16, 32x32, 48x48, 64x64, 128x128, 256x256")
    
    # Clean up temp file
    import os
    os.remove(png_path)
    
except Exception as e:
    print(f"Warning: Could not create icon: {e}")
    print("The executable will use default icon")
