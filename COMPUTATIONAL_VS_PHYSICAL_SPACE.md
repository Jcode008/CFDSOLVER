# Computational vs Physical Space in CFD Solver

## Quick Reference Guide

### **COMPUTATIONAL SPACE** 🔢
Things that happen in **array/index land**:
- Loop indices: `i`, `j` (row, column)
- Grid dimensions: `nx = 120`, `ny = 60`
- Array access: `field.u(i, j)`, `u_star[idx(i,j)]`
- Neighbor offsets: `i±1`, `j±1`
- Index conversion: `idx(i,j) = i*nx + j`

**Think:** "Where am I in the array?"

---

### **PHYSICAL SPACE** 🌍
Things that happen in **real-world units**:
- Velocities: `u`, `v` in **m/s**
- Pressure: `p` in **Pa** (Pascals)
- Cell spacing: `dx`, `dy` in **meters**
- Derivatives: `du/dx` in **1/s**
- Time step: `dt` in **seconds**
- Density: `ρ = 1.225 kg/m³`
- Viscosity: `ν = 1.5×10⁻⁵ m²/s`

**Think:** "What is the actual physical value?"

---

## 🔄 TRANSFORMATION POINTS (Where the Magic Happens!)

These are the **critical lines** where computational differences become physical derivatives:

### 1️⃣ **First Derivative** (Velocity gradient)
```cpp
// COMPUTATIONAL: Access neighbors with index offset
double u_right = field.u(i, j+1);  // Neighbor at [i][j+1]
double u_ij = field.u(i, j);       // Current cell [i][j]

// TRANSFORMATION: Divide difference by PHYSICAL spacing
double du_dx = (u_right - u_ij) / dx;  
//             └─ COMPUTATIONAL ─┘   └─ PHYSICAL ─┘
//                 difference           spacing
// Result: PHYSICAL derivative [m/s]/[m] = [1/s]
```

### 2️⃣ **Second Derivative** (Laplacian/diffusion)
```cpp
// COMPUTATIONAL: 3-point stencil
double u_right = field.u(i, j+1);   // [i][j+1]
double u_ij = field.u(i, j);        // [i][j]
double u_left = field.u(i, j-1);    // [i][j-1]

// TRANSFORMATION: Divide by PHYSICAL spacing squared
double d2u_dx2 = (u_right - 2.0*u_ij + u_left) / (dx * dx);
//               └────── COMPUTATIONAL ──────┘   └─ PHYSICAL ─┘
// Result: PHYSICAL second derivative [m/s]/[m²] = [1/(m·s)]
```

### 3️⃣ **Divergence** (Continuity check)
```cpp
// COMPUTATIONAL: Access neighbors
double u_right = u_star[idx(i, j+1)];
double u_left = u_star[idx(i, j-1)];

// TRANSFORMATION: Central difference divided by spacing
double du_dx = (u_right - u_left) / (2.0 * dx);
// Result: PHYSICAL divergence [1/s]
```

---

## 📊 Complete Data Flow Example

Let's trace one velocity update at cell `[10][20]`:

### **Step 1: Access Data** (Computational → Physical)
```cpp
int i = 10, j = 20;  // COMPUTATIONAL: Array indices

// Access PHYSICAL values stored at these COMPUTATIONAL locations
double u_ij = field.u(10, 20);      // = 5.2 m/s (PHYSICAL)
double u_right = field.u(10, 21);   // = 5.5 m/s (PHYSICAL)
double u_left = field.u(10, 19);    // = 4.9 m/s (PHYSICAL)
```

### **Step 2: Compute Derivative** (TRANSFORMATION!)
```cpp
double dx = 0.0083;  // PHYSICAL spacing [meters]

// Transform COMPUTATIONAL difference → PHYSICAL derivative
double du_dx = (u_right - u_left) / (2.0 * dx);
//           = (5.5 - 4.9) / (2.0 * 0.0083)
//           = 0.6 / 0.0166
//           = 36.1 [1/s]  ← PHYSICAL derivative
```

### **Step 3: Time Integration** (Physical calculation)
```cpp
double dt = 0.001;  // PHYSICAL: [seconds]
double nu = 1.5e-5; // PHYSICAL: [m²/s]

// All PHYSICAL quantities
u_new = u_ij + dt * (-u_ij * du_dx + nu * d2u_dx2);
//      [m/s] + [s] * ([m/s]·[1/s] + [m²/s]·[1/(m·s)])
//      = 5.2 + 0.001 * (-5.2*36.1 + viscous_term)
//      = 5.01 m/s  ← New PHYSICAL velocity
```

### **Step 4: Store Result** (Physical → Computational)
```cpp
// Store PHYSICAL value at COMPUTATIONAL location
field.u(10, 20) = u_new;  // 5.01 m/s stored at [10][20]
```

---

## 🎯 Key Insights

### 1. **Indices are COMPUTATIONAL, Values are PHYSICAL**
```cpp
for (int i = 0; i < ny; i++) {        // i, ny: COMPUTATIONAL
    for (int j = 0; j < nx; j++) {    // j, nx: COMPUTATIONAL
        double u = field.u(i, j);     // u: PHYSICAL [m/s]
        double p = field.p(i, j);     // p: PHYSICAL [Pa]
    }
}
```

### 2. **Spacing Converts Differences to Derivatives**
```cpp
// WRONG: This is just a number difference (dimensionless)
double diff = u[j+1] - u[j];  // Just a difference ❌

// CORRECT: Divide by spacing to get PHYSICAL derivative
double derivative = (u[j+1] - u[j]) / dx;  // [1/s] ✓
//                                    └─ PHYSICAL spacing bridges the gap!
```

### 3. **All Physics Happens in Physical Space**
Every physical equation (Navier-Stokes, continuity, etc.) operates on PHYSICAL quantities:
- ✅ Velocities in m/s
- ✅ Pressure in Pa
- ✅ Derivatives in 1/s or Pa/m
- ✅ Time steps in seconds

The COMPUTATIONAL grid is just a convenient way to **organize** these physical values!

---

## 🧩 Analogy: Spreadsheet vs Real Data

| Aspect | Computational Space | Physical Space |
|--------|-------------------|---------------|
| **What it is** | Excel cell address | Value in the cell |
| **Example** | Row 10, Column 20 (`[10][20]`) | 5.2 m/s |
| **How we use it** | Loop through cells | Calculate physics |
| **Units** | Dimensionless integers | m/s, Pa, meters, seconds |
| **Purpose** | Navigation/storage | Actual simulation |

```
Computational:  [i][j]     →  Cell address "B10"
Physical:       field.u(i,j) →  Value in cell: "5.2 m/s"
Transformation: (u[j+1]-u[j])/dx → Compute rate of change
```

---

## 🔍 Find Transformation Points in Code

Look for these patterns - they're where computational meets physical:

```cpp
// Pattern 1: Division by dx or dy
something / dx        ← Creating a PHYSICAL derivative
something / (dx*dx)   ← Creating a PHYSICAL second derivative

// Pattern 2: Array access followed by division
(u[i][j+1] - u[i][j-1]) / (2.0*dx)  ← Central difference

// Pattern 3: Finite difference stencils
(u[i][j+1] - 2*u[i][j] + u[i][j-1]) / (dx*dx)  ← Laplacian
```

**Rule of thumb:** If you see division by `dx`, `dy`, `dx*dx`, or `dy*dy`, you're at a **transformation point**! 🎯

---

## Summary

- **Computational space** = Index bookkeeping (where to read/write in arrays)
- **Physical space** = Real physics (velocities, pressures, forces in SI units)
- **Transformation** = Dividing differences by physical spacing to get derivatives
- **The bridge** = Cell spacing `dx`, `dy` connects the two worlds

The CFD solver uses **computational indices** to organize data, but performs all **physical calculations** using real-world units! 🚀
