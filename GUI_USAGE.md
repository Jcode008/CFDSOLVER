# CFD Solver GUI - Quick Start Guide

## 🚀 Launching the Application

### On Linux:
```bash
./launch_gui.sh
# or
python3 cfd_gui.py
```

### On Windows:
```powershell
python cfd_gui.py
# or double-click
START_GUI.bat
```

## 📋 How to Use

### Step 1: Set Your Parameters (Simulation Setup Tab)

#### Flow Conditions:
- **Freestream Velocity**: Speed of incoming air (default: 65 m/s = ~126 knots)
- **Angle of Attack**: Airfoil tilt angle in degrees (positive = nose up)
- **Density**: Air density (1.225 kg/m³ = sea level)
- **Kinematic Viscosity**: Air viscosity (1.5e-5 m²/s for standard air)

#### Airfoil Geometry (NACA 4-digit):
- **Chord Length**: Size of airfoil (default: 0.6 m)
- **Max Camber**: NACA first digit / 100 (e.g., 0.02 = 2% = NACA 2xxx)
- **Camber Position**: NACA second digit / 10 (e.g., 0.4 = 40% = NACA x4xx)
- **Max Thickness**: NACA last two digits / 100 (e.g., 0.12 = 12% = NACA xx12)

Example: NACA 2412 → m=0.02, p=0.4, t=0.12

#### Simulation Settings:
- **Grid Points**: Resolution (nx × ny). Higher = better quality but slower
  - 800 × 400 = good quality (~30 mins)
  - 1600 × 800 = high quality (~2 hours)
- **Domain Size**: How far from airfoil (Lx × Ly in meters)
- **Number of Timesteps**: How long to simulate (5000 = ~1 hour real-time)
- **Snapshot Interval**: Save data every N steps (50 = saves 100 snapshots)

#### Click "Calculate Reynolds Number" to check your flow regime:
- Re < 10⁵ = laminar
- Re > 10⁶ = turbulent

#### Click "Apply Parameters & Update Code" to save your settings

### Step 2: Run the Simulation (Run & Monitor Tab)

1. **Enable Live Plotting** (checkbox) if you want real-time visualization
2. **Click "▶ Run Simulation"**
   - The solver will build automatically
   - Progress shows in the console
   - If enabled, a live plot window opens showing velocity and pressure fields
3. **Watch the Progress Bar** - it animates while building/running
4. **Monitor Console Output** - shows timestep progress, convergence info
5. **Stop Button** available if you need to abort

### Step 3: View Results (Results & Analysis Tab)

After simulation completes, results are auto-generated:

- **📊 Generate Analysis**: Creates comprehensive plots with:
  - Velocity field
  - Pressure field
  - Streamlines
  - Lift/Drag coefficients
  
- **🎬 Generate Animation**: Creates GIF showing flow evolution

- **📂 Open Results Folder**: Opens `build/Release/` with all data files

- **🖼️ View Analysis Image**: Opens the summary plot

- **🎥 View Animation**: Opens the flow animation GIF

## 📊 Understanding Results

### Console Output:
```
t=1000/5000 (20.00%)
  CL = 0.86 | CD = 0.013 | L/D = 65.1
  Residual: 1.23e-4
```
- **CL**: Lift coefficient (higher = more lift)
- **CD**: Drag coefficient (lower = less drag)
- **L/D**: Lift-to-drag ratio (efficiency metric)
- **Residual**: Convergence indicator (should decrease)

### Output Files (in `build/Release/`):
- `u_XXXX.csv`, `v_XXXX.csv`: Velocity components at timestep XXXX
- `p_XXXX.csv`: Pressure field
- `metrics_XXXX.csv`: Grid quality metrics
- `aero_forces.csv`: Lift/drag history
- `complete_analysis_100kts.png`: Summary visualization
- `flow_animation_100kts.gif`: Animated flow evolution

## 🎯 Typical Use Cases

### High-Lift Testing (e.g., flaps extended):
```
Velocity: 50 m/s
Angle of Attack: +10°
NACA: 4412 (cambered)
→ Expect CL > 1.0
```

### Cruise Efficiency:
```
Velocity: 80 m/s
Angle of Attack: 0°
NACA: 2412 (moderate camber)
→ Expect CL ≈ 0.3, L/D > 50
```

### Stall Investigation:
```
Velocity: 60 m/s
Angle of Attack: +15° to +20°
→ Watch for flow separation in results
```

## ⚡ Performance Tips

- **Quick Test**: 400×200 grid, 1000 timesteps (~5 mins)
- **Standard Run**: 800×400 grid, 5000 timesteps (~30 mins)
- **High Quality**: 1600×800 grid, 10000 timesteps (~4 hours)

## 🔧 Troubleshooting

### Build Fails:
- Make sure CMake is installed: `cmake --version`
- Check that you have a C++ compiler: `g++ --version`
- Try rebuilding: Click "🔨 Rebuild Only"

### No Live Plot:
- Disable and re-enable the checkbox
- Live plot updates every 2 seconds (be patient)

### Simulation Unstable:
- Reduce velocity
- Increase grid resolution
- Decrease timestep (edit `main.cpp` for advanced users)

### Can't Open Results:
- Check that `xdg-open` works: `xdg-open .`
- Manually navigate to `build/Release/` folder

## 🎨 GUI Features

- ✅ Real-time parameter input with validation
- ✅ Progress bar shows simulation status
- ✅ Live plotting of velocity and pressure fields
- ✅ Automatic result generation
- ✅ Console output for debugging
- ✅ Reynolds number calculator
- ✅ One-click access to results
- ✅ Cross-platform (Windows/Linux/macOS)

## 💡 Tips

1. Always click "Apply Parameters" before running
2. Start with default values for first run
3. Save interesting results (they get overwritten)
4. Watch the console for convergence issues
5. Use live plotting to catch problems early

---

**Enjoy your CFD simulations!** 🚁✈️
