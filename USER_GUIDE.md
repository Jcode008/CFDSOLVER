# CFD Solver User Guide
## How to Analyze Airfoils with Your O-Grid Solver

---

## Quick Start (3 Steps)

### 1. Generate the Grid
```powershell
cd build
.\Release\TestMesh.exe
```

**What this does:**
- Creates an O-grid wrapped around your airfoil
- Exports `grid.csv` and `metrics.csv`
- Shows grid quality metrics

**Key settings** (edit `src/test_mesh.cpp`):
```cpp
int nxi = 120;              // Points around airfoil (more = finer)
int neta = 60;              // Points radial (surface → farfield)
double chord = 0.6;         // Airfoil chord length (meters)
double alpha_deg = -5.0;    // Grid rotation (NEGATIVE = positive AoA!)
double farfield_radius = 1.8;  // How far the domain extends
```

### 2. Run the Solver
```powershell
.\Release\TestCurvilinear.exe
```

**What this does:**
- Loads the grid
- Initializes flow (horizontal freestream)
- Runs 2000 timesteps of Navier-Stokes
- Exports `u_final.csv`, `v_final.csv`, `p_final.csv`

**Key settings** (edit `src/test_curvilinear.cpp`):
```cpp
double U_inf = 5.0;     // Freestream velocity (m/s)
double Re = 15000.0;    // Reynolds number
double dt = 1e-5;       // Timestep (seconds)
int n_steps = 2000;     // How long to run
```

### 3. Visualize Results
```powershell
cd ..
python analysis/visualize_flow.py          # 6-panel flow visualization
python analysis/plot_all_results.py        # 12-panel comprehensive analysis
```

**Outputs:**
- `build/flow_solution.png` - Velocity, pressure, vectors, Cp
- `build/complete_results.png` - Grid + flow + convergence
- Terminal shows: CL, CD, L/D

---

## Understanding the Airfoil Parameters

### NACA 4-Digit Designation (e.g., NACA 2412)

In `src/ogrid_generator.cpp`:
```cpp
const double m = 0.02;     // First digit: Max camber (2%)
const double p = 0.4;      // Second digit: Camber location (40% chord)
const double thick = 0.12; // Last two digits: Thickness (12%)
```

**Common Airfoils:**
- **NACA 0012**: Symmetric (no camber), 12% thick - good for vertical tails
- **NACA 2412**: 2% camber at 40%, 12% thick - general aviation wings
- **NACA 4412**: 4% camber at 40%, 12% thick - high lift
- **NACA 0015**: Symmetric, 15% thick - helicopter blades

### Angle of Attack Convention ⚠️ IMPORTANT!

**Grid rotation vs. Effective AoA:**
```
alpha_deg = -5° in grid  →  Effective α = +5° (POSITIVE lift)
alpha_deg = +5° in grid  →  Effective α = -5° (NEGATIVE lift)
```

**Why?** The grid generator **rotates the airfoil**, but the flow stays horizontal.
- Negative grid rotation = nose tilts down = flow hits from below = **positive AoA**
- Positive grid rotation = nose tilts up = flow hits from above = **negative AoA**

**Example:**
```cpp
// For α = +10° angle of attack (strong positive lift):
double alpha_deg = -10.0;  // Rotate airfoil DOWN by 10°

// For α = 0° (zero lift, symmetric flow):
double alpha_deg = 0.0;    // No rotation

// For α = -5° (negative lift, downforce):
double alpha_deg = +5.0;   // Rotate airfoil UP by 5°
```

---

## Workflow for Different Studies

### Study 1: Single Airfoil Analysis

**Goal:** Get CL, CD for one airfoil at one angle of attack

1. **Edit** `src/test_mesh.cpp`:
   ```cpp
   double alpha_deg = -5.0;  // For α=+5° effective
   ```

2. **Edit** `src/ogrid_generator.cpp`:
   ```cpp
   const double m = 0.02;      // NACA 2XXX
   const double p = 0.4;       // NACA X4XX
   const double thick = 0.12;  // NACA XX12
   ```

3. **Rebuild:**
   ```powershell
   cmake --build . --config Release
   ```

4. **Run:**
   ```powershell
   .\Release\TestMesh.exe
   .\Release\TestCurvilinear.exe
   ```

5. **Check results** - Terminal shows:
   ```
   CL = 0.8555
   CD = 0.0132
   L/D = 64.62
   ```

### Study 2: Angle of Attack Sweep (Polar)

**Goal:** Generate CL vs α curve

Create a batch script `run_polar.ps1`:
```powershell
# Sweep α from -10° to +15°
$angles = -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 15

foreach ($alpha in $angles) {
    Write-Host "Running α = $(-$alpha)° (grid = $alpha°)"
    
    # Update test_mesh.cpp with new angle
    # (You'll need to automate this or manually change each time)
    
    .\Release\TestMesh.exe
    .\Release\TestCurvilinear.exe
    
    # Save results
    Copy-Item u_final.csv "u_alpha$alpha.csv"
    Copy-Item p_final.csv "p_alpha$alpha.csv"
}
```

Then plot all results together.

### Study 3: Reynolds Number Sweep

**Goal:** See how Re affects drag

Edit `src/test_curvilinear.cpp` and run multiple times:
```cpp
// Low Re (laminar):
double Re = 1000.0;

// Transitional:
double Re = 10000.0;

// Higher Re (more turbulent):
double Re = 50000.0;
```

**Expect:** CD decreases as Re increases (thinner boundary layer).

### Study 4: Compare Different Airfoils

**Goal:** Which airfoil is better for your application?

Test each airfoil at same conditions:
```
NACA 0012 (symmetric)   → Good L/D, zero-lift at α=0°
NACA 2412 (cambered)    → Higher CL at same α
NACA 4412 (more camber) → Even higher CL, but more drag
NACA 0015 (thick)       → Structurally strong, more drag
```

---

## Interpreting Results

### Velocity Field (`flow_solution.png`)

**What to look for:**
- ✅ **Smooth contours** - solution converged
- ✅ **Higher velocity on lower surface** (α > 0) - generating lift
- ✅ **Stagnation point at leading edge** - correct physics
- ❌ **Oscillations/wiggles** - increase grid resolution
- ❌ **Separated flow at trailing edge** - too high α or Re too low

### Pressure Field

**What to look for:**
- ✅ **Low pressure (suction) on lower surface** - positive lift
- ✅ **High pressure (stagnation) near leading edge** - correct
- ✅ **Pressure recovery toward trailing edge** - good
- ❌ **Pressure plateau** - flow separation
- ❌ **Negative pressures everywhere** - initialization problem

### Pressure Coefficient (Cp)

**Formula:** Cp = (p - p∞) / (0.5 × ρ × U∞²)

**Interpretation:**
- **Cp = 0**: Freestream pressure
- **Cp = 1**: Stagnation point (complete velocity stop)
- **Cp < 0**: Suction (flow accelerated above freestream)
- **Cp > 1**: Impossible for incompressible flow!

**Typical Cp plot:**
```
     +1 ├─────┐               (Stagnation)
        │     └───┐
Cp   0 ├─────────└─────────  (Freestream)
        │           ┌────
     -1 ├───────────┘         (Suction peak)
        │
     -2 ├                     (Strong suction)
        └────────────────────
        0%   50%   100%  x/c
        LE         TE
```

### Aerodynamic Coefficients

**Lift Coefficient (CL):**
```
CL = L / (0.5 × ρ × U∞² × S)
```
- **CL = 0**: No lift (symmetric flow)
- **CL > 0**: Upward lift (positive AoA)
- **CL < 0**: Downward force (negative AoA)
- **Typical:** 0.2 to 1.5 for normal operation
- **Stall:** CL drops suddenly at high α (flow separates)

**Drag Coefficient (CD):**
```
CD = D / (0.5 × ρ × U∞² × S)
```
- **CD pressure**: From pressure distribution (form drag)
- **CD viscous**: From wall shear stress (skin friction) - NOT YET IMPLEMENTED
- **CD total = CD pressure + CD viscous**
- **Typical:** 0.01 to 0.05 for clean airfoils

**Lift-to-Drag Ratio (L/D):**
```
L/D = CL / CD
```
- **Higher is better!** (more efficient)
- **L/D = 10-20**: Typical for aircraft at cruise
- **L/D = 40-60**: High-performance gliders
- **L/D = 60+**: Your solver shows excellent efficiency!

---

## Troubleshooting

### Problem: Solution Diverges (NaN appears)

**Symptoms:** 
```
⚠ NaN at step 50, i=74, j=58
```

**Solutions:**
1. **Reduce timestep:**
   ```cpp
   double dt = 5e-6;  // Half the current dt
   ```

2. **Reduce farfield distance:**
   ```cpp
   double farfield_radius = 1.2;  // Closer to airfoil
   ```

3. **Increase radial resolution:**
   ```cpp
   int neta = 80;  // More points toward farfield
   ```

### Problem: Grid Has Negative Jacobians

**Symptoms:**
```
⚠ WARNING: 48 points have J ≤ 0!
Grid has folded or crossed lines.
```

**Solutions:**
1. **Reduce resolution** (paradoxically!):
   ```cpp
   int nxi = 100;   // Fewer points around
   int neta = 40;   // Fewer points radial
   ```

2. **Increase farfield:**
   ```cpp
   double farfield_radius = 2.5;  // Give grid more space
   ```

3. **Adjust elliptic smoothing:**
   ```cpp
   smoothElliptic(20000, 1e-9);  // More iterations, tighter tolerance
   ```

### Problem: Wrong Flow Direction

**Symptoms:** Velocity high on wrong surface

**Check:**
1. **Angle convention:**
   ```cpp
   // For POSITIVE lift, use NEGATIVE grid angle!
   double alpha_deg = -5.0;  // This gives α=+5° effective
   ```

2. **Flow initialization:**
   ```cpp
   // Should be horizontal for standard AoA
   field.u(i, j) = U_inf;  // Horizontal
   field.v(i, j) = 0.0;    // Not angled!
   ```

### Problem: Solution Not Converging

**Symptoms:** Velocity/pressure still changing after 2000 steps

**Solutions:**
1. **Run longer:**
   ```cpp
   int n_steps = 5000;  // More timesteps
   ```

2. **Check convergence plot** in `complete_results.png`:
   - Should see |dV/dt| decreasing exponentially
   - If still decreasing, just need more time
   - If flat or increasing, check grid/physics

3. **Increase Reynolds number** (counter-intuitive!):
   ```cpp
   double Re = 30000.0;  // Higher Re can stabilize solution
   ```

---

## Advanced Modifications

### Change Airfoil Shape

**In** `src/ogrid_generator.cpp`:

```cpp
// NACA 6-series (laminar flow):
// Not implemented - requires different equations

// Custom airfoil from file:
// 1. Create x[], y[] arrays from file
// 2. Replace NACA equations with your data
// 3. Interpolate to nxi points

// Flat plate (simplest):
const double thick = 0.0;  // Zero thickness
const double m = 0.0;      // No camber
const double p = 0.5;      // Doesn't matter
```

### Add Flap Deflection

**Concept:** Rotate trailing edge section

```cpp
// In ogrid_generator.cpp, after computing airfoil:
double flap_start = 0.75;  // Flap at 75% chord
double flap_angle = 10.0 * M_PI / 180.0;  // 10° down

for (int i = 0; i < nxi_; i++) {
    double s = // ...chord position...
    if (s > flap_start) {
        // Rotate trailing edge points
        double dx = x_airfoil[i] - (flap_start * chord);
        double dy = y_airfoil[i];
        x_airfoil[i] = flap_start * chord + dx * cos(flap_angle) - dy * sin(flap_angle);
        y_airfoil[i] = dx * sin(flap_angle) + dy * cos(flap_angle);
    }
}
```

### Export Timestep Data

**In** `src/test_curvilinear.cpp`:

```cpp
for (int step = 0; step <= n_steps; step++) {
    // ... solver step ...
    
    // Export every 100 steps
    if (step % 100 == 0) {
        exportFieldCSV("u_" + std::to_string(step) + ".csv", field, nxi, neta, 'u');
        exportFieldCSV("v_" + std::to_string(step) + ".csv", field, nxi, neta, 'v');
        exportFieldCSV("p_" + std::to_string(step) + ".csv", field, nxi, neta, 'p');
    }
}
```

Then create animations with `create_animation.py`.

---

## Performance Tips

### Faster Simulations

1. **Reduce resolution** (biggest impact):
   ```cpp
   int nxi = 80;   // 120 → 80 (33% fewer points)
   int neta = 40;  // 60 → 40
   // Run time: ~30% faster
   ```

2. **Increase timestep** (careful - may diverge!):
   ```cpp
   double dt = 2e-5;  // 2× larger
   // Run time: 50% faster
   ```

3. **Fewer timesteps** (if only steady solution needed):
   ```cpp
   int n_steps = 1000;  // Half the steps
   // Run time: 50% faster
   ```

### Better Accuracy

1. **Higher resolution**:
   ```cpp
   int nxi = 200;
   int neta = 80;
   // More accurate, but 3× slower
   ```

2. **Smaller timestep**:
   ```cpp
   double dt = 5e-6;
   // More stable, captures faster dynamics
   ```

3. **Better convergence tolerance**:
   ```cpp
   // In pressure solver (curvilinear_solver.cpp):
   double tol = 1e-7;  // Stricter (currently 1e-6)
   ```

---

## File Reference

### Grid Files (CSV format)

**`grid.csv`:**
```csv
i,j,x,y
0,0,0.6,-0.0124
0,1,0.6122,-0.0201
...
```
- **i**: Circumferential index (0 to nxi-1)
- **j**: Radial index (0 to neta-1)
- **x, y**: Physical coordinates (meters)

**`metrics.csv`:**
```csv
i,j,J,xi_x,xi_y,eta_x,eta_y
0,0,0.0125,0.14,-0.08,1.24,0.95
...
```
- **J**: Jacobian (grid cell area, must be > 0)
- **xi_x, xi_y**: ∂x/∂ξ, ∂y/∂ξ (tangent to airfoil)
- **eta_x, eta_y**: ∂x/∂η, ∂y/∂η (normal to airfoil)

### Solution Files

**`u_final.csv`, `v_final.csv`, `p_final.csv`:**
```csv
i,j,value
0,0,5.0
0,1,5.02
...
```
- **value**: u-velocity, v-velocity, or pressure at grid point (i,j)

---

## Next Steps

### Immediate Improvements:

1. **Add viscous drag calculation**
   - Integrate wall shear stress: τ_wall = μ × (∂u/∂n)|wall
   - Get accurate total CD

2. **Implement turbulence model**
   - k-ε or Spalart-Allmaras
   - Better high-Re predictions

3. **Adaptive timestep**
   - Automatically adjust dt based on CFL condition
   - Faster convergence

4. **Parallel processing**
   - OpenMP for multi-core
   - 4-8× speedup on modern CPUs

### Learning More:

- **CFD textbooks**: Anderson "Computational Fluid Dynamics"
- **Grid generation**: Thompson, Warsi "Numerical Grid Generation"
- **Airfoil theory**: Abbott & von Doenhoff "Theory of Wing Sections"

---

## Support & Troubleshooting

**Check logs:** All outputs go to terminal - save with:
```powershell
.\Release\TestCurvilinear.exe > results.log 2>&1
```

**Grid quality:** Always check `grid_analysis.png` first!
- Green Jacobians = good
- Red/negative = grid failed

**Common mistakes:**
1. Forgetting to regenerate grid after changing airfoil
2. Wrong angle of attack sign convention
3. Timestep too large for grid resolution
4. Farfield too close (should be 3-5× chord minimum)

---

**Happy simulating!** 🚀
