# CFD Solver GUI - Standalone Executable

## 🚀 Quick Start

### Option 1: Build the Executable (Recommended)

Simply run:
```bash
build_exe.bat
```

This will:
- Check if PyInstaller is installed (installs it if needed)
- Build a standalone `CFDSolverGUI.exe` in the `dist` folder
- Takes 2-5 minutes to complete

### Option 2: Run with Python (Development)

```bash
python cfd_gui.py
```

## 📦 What You Get

After building, you'll have:
- `dist/CFDSolverGUI.exe` - **Standalone executable** (~150-200 MB)
  - No Python installation required to run!
  - Bundles Python + NumPy + Matplotlib + all dependencies
  - Can be shared with others

## ⚙️ System Requirements

### To Run the .exe:
- Windows 10/11 (64-bit)
- ~4 GB RAM minimum (8+ GB recommended)
- ~2 GB free disk space for results
- CMake and C++ build tools (for the solver itself)

### To Build the .exe:
- Python 3.8+ with packages:
  - tkinter (usually included)
  - numpy
  - matplotlib
  - pyinstaller

Install dependencies:
```bash
pip install numpy matplotlib pyinstaller
```

## 🎯 Features

### GUI Application:
✅ **Interactive Parameter Setup**
- Flow conditions (velocity, angle of attack, viscosity)
- Airfoil geometry (NACA 4-digit parameters)
- Grid resolution and simulation settings
- Reynolds number calculator

✅ **One-Click Execution**
- Automatic build and run
- Real-time console output
- Progress monitoring
- Stop/restart capability

✅ **Live Flow Visualization** 🔥
- Real-time velocity and pressure plots
- Updates every 2 seconds during simulation
- Separate window with dual plots
- Timestep and statistics display

✅ **Automatic Post-Processing**
- Animation generation (GIF)
- Comprehensive analysis plots
- Quick-view buttons
- Results folder access

## 📊 Using the Application

### 1. Launch
```bash
dist\CFDSolverGUI.exe
```
or double-click `CFDSolverGUI.exe`

### 2. Setup Simulation (Tab 1)
- Enter your parameters
- Click "Calculate Reynolds Number"
- Click "Apply Parameters"

### 3. Run (Tab 2)
- Check "Enable Live Plotting" (recommended!)
- Click "▶ Run Simulation"
- Watch live plots update in real-time
- Monitor console for progress

### 4. View Results (Tab 3)
- Animation and analysis auto-generate
- Click buttons to view results
- Open results folder for CSV data

## 🎨 Live Plotting

When enabled (checkbox in Tab 2), a separate window opens showing:
- **Left plot**: Velocity magnitude field (color-coded)
- **Right plot**: Pressure field (red/blue diverging)
- **Status bar**: Current timestep, physical time, min/max values
- **Auto-updates**: Every 2 seconds with latest snapshot

## 📁 File Structure

```
CFDSolver/
├── cfd_gui.py              # Main GUI application
├── build_exe.bat           # One-click executable builder
├── setup.py                # Alternative PyInstaller setup
├── dist/
│   └── CFDSolverGUI.exe   # Standalone executable (after build)
├── build/
│   └── Release/
│       ├── CFDSolver.exe  # C++ solver
│       ├── u_*.csv        # Velocity data
│       ├── p_*.csv        # Pressure data
│       └── results...     # Animations, plots
├── analysis/              # Python post-processing scripts
├── src/                   # C++ source code
└── include/              # C++ headers
```

## 🔧 Troubleshooting

### "PyInstaller not found"
```bash
pip install pyinstaller
```

### "Module not found" when running .exe
Rebuild with:
```bash
build_exe.bat
```

### Live plotting not updating
- Check that simulation is running
- Ensure snapshot files are being created in `build/Release/`
- Verify checkbox is enabled before clicking Run

### Simulation won't start
- Make sure CMake is configured: `cmake -B build`
- Check that C++ solver builds: `cmake --build build --config Release`

## 📝 Distribution

To share your application:

1. **Copy these files:**
   - `dist/CFDSolverGUI.exe`
   - `build/` folder (contains C++ solver)
   - `analysis/` folder (post-processing scripts)

2. **Or create an installer:**
   - Use Inno Setup or NSIS
   - Include CMake dependencies
   - Package everything together

3. **GitHub Release:**
   - Upload `CFDSolverGUI.exe` as a release asset
   - Include this README
   - Tag with version number

## 🌟 Publishing Checklist

- [ ] Test executable on clean Windows machine
- [ ] Create GitHub repository with README
- [ ] Add screenshots of GUI and results
- [ ] Write detailed usage documentation
- [ ] Create example parameter sets
- [ ] Add license (MIT recommended)
- [ ] Create release with .exe file
- [ ] Consider creating video tutorial
- [ ] Share on relevant forums/communities

## 📧 Support

For issues or questions:
- Check console output for errors
- Verify system requirements
- Review this README thoroughly
- Open an issue on GitHub

## 🎓 Citation

If you use this software in research or publications, please cite:
```
CFD Solver GUI - Interactive Airfoil Flow Simulation
[Your Name/Organization]
[Year]
[GitHub URL]
```

## 📄 License

[Choose your license - MIT recommended for open source]

---

**Enjoy your CFD simulations!** 🚀
