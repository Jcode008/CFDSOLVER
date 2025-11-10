# CFD Solver for Airfoil Analysis

A body-fitted computational fluid dynamics (CFD) solver for analyzing 2D airfoils using O-grid topology and curvilinear coordinates.

## 🚀 Features

- **O-Grid Generation**: Body-fitted mesh wrapped around airfoils (no wake cut!)
- **NACA 4-Digit Airfoils**: Built-in NACA airfoil geometry generator
- **Fully Implicit Navier-Stokes**: Crank-Nicolson diffusion + upwind advection
- **Stable Long-Duration Runs**: 2000+ timesteps (40× improvement over C-grid)
- **Comprehensive Visualization**: Grid quality, flow field, aerodynamic analysis

## 📊 Results

**NACA 2412 at α=0°:**
- ✅ Lift coefficient: CL ≈ 0.3
- ✅ Drag coefficient: CD ≈ 0.01
- ✅ L/D ratio: ~30
- ✅ Stable for 2000+ timesteps

**NACA 2412 at α=+5°:**
- ✅ Lift coefficient: CL ≈ 0.86
- ✅ Drag coefficient: CD ≈ 0.013
- ✅ L/D ratio: ~65

## 🛠️ Build Instructions

### Prerequisites
- CMake 3.10+
- C++ compiler (MSVC, GCC, or Clang)
- Python 3.x with: `numpy`, `matplotlib`, `pandas`

### Build
```powershell
mkdir build
cd build
cmake ..
cmake --build . --config Release
```

## 🎯 Quick Start

### 1. Generate O-Grid
```powershell
cd build
.\Release\TestMesh.exe
```

### 2. Run CFD Solver
```powershell
.\Release\TestCurvilinear.exe
```

### 3. Visualize Results
```powershell
cd ..
python analysis/visualize_flow.py
python analysis/plot_all_results.py
```

## 📁 Project Structure

```
CFDSolver/
├── src/                    # C++ source files
│   ├── main.cpp           # Basic solver (unused)
│   ├── test_mesh.cpp      # O-grid generator driver
│   ├── test_curvilinear.cpp  # Curvilinear solver driver
│   ├── grid.cpp           # Grid data structures
│   ├── field.cpp          # Flow field
│   ├── solver.cpp         # Navier-Stokes solver
│   └── ogrid_generator.cpp   # O-grid generation
├── include/               # Header files
│   ├── grid.hpp
│   ├── field.hpp
│   ├── solver.hpp
│   └── utils.hpp
├── analysis/              # Python visualization scripts
│   ├── visualize_flow.py
│   ├── visualize_grid.py
│   └── plot_all_results.py
├── build/                 # Build output (gitignored)
├── CMakeLists.txt        # CMake configuration
├── USER_GUIDE.md         # Comprehensive user guide
└── README.md             # This file
```

## 🎓 How It Works

### O-Grid Topology
Unlike traditional C-grids (which have a wake cut), O-grids wrap completely around the airfoil:
- **Smoother grid metrics** - no discontinuity at wake
- **Better aspect ratio** - 116k:1 vs 2.3M:1 for C-grid
- **Stable simulations** - runs 40× longer before divergence

### Curvilinear Coordinate Transform
Flow equations solved in computational space (ξ,η) using contravariant velocity formulation:
```
U = u·ξ_x + v·ξ_y
V = u·η_x + v·η_y
```

### Numerical Scheme
- **Time integration**: Fractional step (projection method)
- **Diffusion**: Crank-Nicolson (implicit, unconditionally stable)
- **Advection**: Upwind scheme (implicit for stability)
- **Pressure**: Poisson equation solved with Gauss-Seidel

### Grid Generation
1. **Algebraic initialization**: Radial lines from airfoil to farfield
2. **Elliptic smoothing**: Laplace equation for smooth distribution
3. **Metric calculation**: Jacobian, transformation derivatives

## 🔧 Customization

### Change Airfoil
Edit `src/ogrid_generator.cpp`:
```cpp
const double m = 0.02;     // Max camber (2% for NACA 2412)
const double p = 0.4;      // Camber location (40% chord)
const double thick = 0.12; // Thickness (12%)
```

### Change Angle of Attack
Edit `src/test_mesh.cpp`:
```cpp
double alpha_deg = -5.0;  // Grid rotation
// NOTE: α_grid = -5° gives effective α = +5° (opposite sign!)
```

### Adjust Grid Resolution
Edit `src/test_mesh.cpp`:
```cpp
int nxi = 120;   // Points around airfoil
int neta = 60;   // Points radial (surface → farfield)
```

### Solver Parameters
Edit `src/test_curvilinear.cpp`:
```cpp
double U_inf = 5.0;    // Freestream velocity (m/s)
double Re = 15000.0;   // Reynolds number
double dt = 1e-5;      // Timestep (s)
int n_steps = 2000;    // Number of timesteps
```

## 📖 Documentation

See **[USER_GUIDE.md](USER_GUIDE.md)** for:
- Detailed workflow tutorials
- Result interpretation guide
- Troubleshooting common issues
- Advanced modifications
- Performance tips

## 🐛 Troubleshooting

**Solution diverges (NaN)?**
- Reduce timestep: `dt = 5e-6`
- Increase radial resolution: `neta = 80`

**Grid has negative Jacobians?**
- Reduce resolution (paradoxically helps!)
- Increase farfield distance

**Wrong flow direction?**
- Check angle sign: grid rotation is OPPOSITE of effective AoA
- For positive lift: use NEGATIVE grid angle

## 📚 References

- Anderson, J.D. - *Computational Fluid Dynamics: The Basics with Applications*
- Thompson, Warsi, Mastin - *Numerical Grid Generation*
- Abbott & von Doenhoff - *Theory of Wing Sections*

## 🤝 Contributing

This is a research/educational project. Contributions welcome!

## 📄 License

MIT License - feel free to use for learning and research.

## ✨ Acknowledgments

Built with blood, sweat, and a lot of debugging. Special thanks to whoever invented implicit schemes.

---

**Want to learn more?** Check out `USER_GUIDE.md` for the complete walkthrough!
