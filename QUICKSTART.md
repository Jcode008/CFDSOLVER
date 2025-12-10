# 🚀 Quick Start Guide - CFD Solver GUI

## What's New? 🎉

### ✨ Live Plotting Feature
- **Real-time visualization** of flow evolution
- Opens in separate window during simulation
- Shows velocity magnitude and pressure fields
- Updates every 2 seconds with latest snapshot
- Displays current timestep, physical time, and statistics

### 📦 Standalone Executable
- Create `CFDSolverGUI.exe` - runs without Python!
- Share with anyone (no installation required)
- One-click build process

---

## How to Use Live Plotting

1. **Launch the GUI** (the window should be open now)
2. **Go to Tab 2** (Run & Monitor)
3. **Check the box** "📈 Enable Live Plotting"
4. **Click "▶ Run Simulation"**
5. **Watch!** A new window opens showing live flow fields

### What You'll See:
```
┌─────────────────────────────────────┐
│  Live Flow Visualization            │
├──────────────────┬──────────────────┤
│  Velocity        │  Pressure        │
│  Magnitude       │  Field           │
│  (colorful!)     │  (red/blue)      │
│                  │                  │
│  [flow plot]     │  [pressure plot] │
│                  │                  │
└──────────────────┴──────────────────┘
│ Timestep: 150 | Time: 0.0019s | ... │
└──────────────────────────────────────┘
```

**Updates every 2 seconds** with newest data!

---

## How to Build Standalone .exe

### Method 1: Super Easy (Recommended) 🎯

Just double-click:
```
build_exe.bat
```

That's it! Wait 2-5 minutes. You'll get `dist/CFDSolverGUI.exe`

### Method 2: Command Line

```bash
cd C:\Users\graha\CFDSolver
build_exe.bat
```

Or manually:
```bash
pyinstaller cfd_gui.py --name=CFDSolverGUI --onefile --windowed
```

### What Happens During Build:

1. ✓ Checks for PyInstaller (installs if needed)
2. ✓ Analyzes Python script and dependencies
3. ✓ Bundles Python interpreter + all libraries
4. ✓ Packages everything into single .exe
5. ✓ Creates `dist/CFDSolverGUI.exe` (~150-200 MB)

**Build time:** 2-5 minutes
**Result:** Fully standalone application!

---

## Testing the .exe

After building:

```bash
cd dist
CFDSolverGUI.exe
```

Or just double-click `CFDSolverGUI.exe` in File Explorer!

**No Python needed to run the .exe!** 🎉

---

## Distribution Tips

### For Friends/Colleagues:
1. Give them `CFDSolverGUI.exe`
2. They also need the `build/` folder (contains C++ solver)
3. And the `analysis/` folder (post-processing scripts)

### For GitHub Release:
1. Create a release tag (e.g., `v1.0.0`)
2. Upload `CFDSolverGUI.exe` as asset
3. Include README with instructions
4. Add screenshots/demo GIF

### For Publication:
1. Create installer with Inno Setup
2. Package .exe + solver + examples
3. Add documentation PDF
4. Include sample datasets

---

## Current Simulation Parameters

Check the GUI - Tab 1 shows all current values:
- Velocity: 65 m/s
- Angle of attack: 0°
- Viscosity: 1.0e-5 m²/s
- NACA 2412 airfoil
- Grid: 800 × 400
- 5000 timesteps

**Try changing them and running!**

---

## Pro Tips 💡

### For Faster Simulations:
- Reduce grid: 400 × 200
- Fewer timesteps: 1000
- Larger snapshot interval: 100

### For Better Quality:
- Increase grid: 1000 × 500
- More timesteps: 10000
- Smaller snapshot interval: 25

### For Cool Visualizations:
- Enable live plotting ✓
- Try different angles of attack (0°, 5°, 10°)
- Vary Reynolds number (change velocity or viscosity)

---

## Troubleshooting

### "Live plot not opening"
- Make sure checkbox is checked **before** clicking Run
- Check that simulation is actually running (console output)

### "Exe build failed"
- Install missing packages: `pip install numpy matplotlib`
- Try running `build_exe.bat` again

### "Simulation crashes"
- Check parameters are reasonable
- Ensure CMake build succeeded
- Look at console output for errors

---

## Next Steps

1. **Test the live plotting** - Run a simulation now!
2. **Build the .exe** - Run `build_exe.bat`
3. **Try different parameters** - Play with angles of attack
4. **Share your results** - Post animations/plots
5. **Publish!** - GitHub, paper, presentation

**Enjoy!** 🎊
