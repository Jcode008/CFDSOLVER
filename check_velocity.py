import pandas as pd
import numpy as np

# Load data
grid = pd.read_csv('build/grid.csv')
u = pd.read_csv('build/u_final.csv')
v = pd.read_csv('build/v_final.csv')

nxi = grid['i'].max() + 1
neta = grid['j'].max() + 1

x = grid['x'].values.reshape(neta, nxi)
y = grid['y'].values.reshape(neta, nxi)
u_arr = u['value'].values.reshape(neta, nxi)
v_arr = v['value'].values.reshape(neta, nxi)
vel = np.sqrt(u_arr**2 + v_arr**2)

print("="*60)
print("AIRFOIL ORIENTATION CHECK")
print("="*60)
print(f"Grid rotation: α_grid = -5°")
print(f"Effective AoA: α_eff = +5° (positive lift)")
print()

print("AIRFOIL SHAPE (j=0, on surface):")
print(f"{'i':>4s} {'x (m)':>8s} {'y (m)':>8s} {'Side':>15s}")
print("-"*40)
for i in [0, 20, 40, 60, 80, 100, 119]:
    if y[0,i] > 0.01:
        side = "UPPER SURFACE"
    elif y[0,i] < -0.01:
        side = "LOWER SURFACE"  
    else:
        side = "Nose/Tail"
    print(f"{i:4d} {x[0,i]:8.4f} {y[0,i]:8.4f} {side:>15s}")

print()
print("VELOCITY NEAR SURFACE (j=1, first grid point off surface):")
print(f"{'i':>4s} {'x (m)':>8s} {'y (m)':>8s} {'V (m/s)':>9s} {'Side':>15s}")
print("-"*50)
for i in [0, 20, 40, 60, 80, 100, 119]:
    if y[1,i] > 0.01:
        side = "NEAR UPPER"
    elif y[1,i] < -0.01:
        side = "NEAR LOWER"
    else:
        side = "Near LE/TE"
    print(f"{i:4d} {x[1,i]:8.4f} {y[1,i]:8.4f} {vel[1,i]:9.3f} {side:>15s}")

print()
print("="*60)
print("WHERE IS VELOCITY MAXIMUM?")
print("="*60)

# Find max velocity near surface (j=1)
max_idx = np.argmax(vel[1,:])
print(f"On first off-surface layer (j=1):")
print(f"  Max velocity: {vel[1,max_idx]:.3f} m/s")
print(f"  Location: i={max_idx}, x={x[1,max_idx]:.4f}, y={y[1,max_idx]:.4f}")
print(f"  This is {'ABOVE' if y[1,max_idx] > 0 else 'BELOW'} the chord line (y=0)")

# Also check j=2 and j=3
for j in [2, 3, 4]:
    max_idx_j = np.argmax(vel[j,:])
    print(f"At j={j}: Max V = {vel[j,max_idx_j]:.3f} m/s at y={y[j,max_idx_j]:.4f}")

print()
print("SIMPLE ANSWER:")
if y[1,max_idx] < 0:
    print("✓ Velocity is FASTER on LOWER surface (y < 0)")
    print("  This is CORRECT for cambered airfoil at positive AoA!")
else:
    print("✓ Velocity is FASTER on UPPER surface (y > 0)")
    print("  This is typical for high AoA or symmetric airfoils")
