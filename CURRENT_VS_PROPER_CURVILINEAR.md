# Current Code vs Proper Curvilinear Implementation

## 🚨 CRITICAL UNDERSTANDING 🚨

### What the Current Code Does (CARTESIAN GRID)

The current `solver.cpp` is treating the O-grid **AS IF it were a simple rectangular Cartesian grid**!

```cpp
// Current code does this (WRONG for O-grid):
double du_dx = (u_right - u_left) / (2.0 * dx);  // ❌ Assumes dx is constant everywhere
```

**Problem:** 
- We use `dx = 0.0083 m` (constant) everywhere
- But in the actual O-grid, cells are NOT uniform rectangles!
- Near the airfoil, cells are tiny and stretched
- Far from airfoil, cells are larger
- **We're pretending the curved grid is a straight rectangle!**

---

## What SHOULD Happen (Proper Curvilinear Grid Transformation)

### The Two Approaches:

## **Approach 1: Transform Derivatives (What you described)** ✓

**YES, you're exactly right!** This is the proper way to handle curvilinear grids:

### Step 1: Work in Computational Space (Rectangle)
```
Computational Space (ξ,η):
┌─────────────────┐  ← Perfect rectangle
│  [i][j] cells   │  ← Uniform spacing: Δξ = 1, Δη = 1
│  Solve here!    │  ← All equations solved on this rectangle
└─────────────────┘
```

**In computational space:**
- Grid is perfectly rectangular (120×60)
- All cells are identical squares
- Spacing is uniform: `Δξ = 1`, `Δη = 1`
- We solve pressure, velocity on THIS rectangle

### Step 2: Compute Derivatives Using Chain Rule + Jacobian

When we need physical derivatives ∂u/∂x, we use the transformation:

```
∂u/∂x = ξ_x · ∂u/∂ξ + η_x · ∂u/∂η

where:
  ξ_x = ∂ξ/∂x  }
  ξ_y = ∂ξ/∂y  } ← These are the METRIC TERMS (inverse Jacobian)
  η_x = ∂η/∂x  }   computed from grid geometry
  η_y = ∂η/∂y  }
```

**Example for ∂u/∂x:**
```cpp
// Step 1: Compute computational derivatives (easy - uniform spacing!)
double du_dxi = (u[i][j+1] - u[i][j-1]) / (2.0 * 1.0);  // Δξ = 1
double du_deta = (u[i+1][j] - u[i-1][j]) / (2.0 * 1.0); // Δη = 1

// Step 2: Transform to physical using metrics
double du_dx = metrics.xi_x[i][j] * du_dxi + metrics.eta_x[i][j] * du_deta;
//             └──────── JACOBIAN MULTIPLIERS ────────┘
```

### Step 3: Solve Equations in Computational Space

**Navier-Stokes in transformed coordinates:**
```cpp
// Continuity (incompressibility): ∇·v = 0 becomes
∂(U)/∂ξ + ∂(V)/∂η = 0

where U, V are contravariant velocities:
  U = J(ξ_x·u + ξ_y·v)  ← Transformed velocity components
  V = J(η_x·u + η_y·v)

// Pressure Poisson becomes:
∂/∂ξ[J²(ξ_x² + ξ_y²)∂p/∂ξ + J²(ξ_x·η_x + ξ_y·η_y)∂p/∂η] + ...
└──────────── Metric terms everywhere! ────────────┘
```

**Key point:** The Jacobian J and metrics (ξ_x, ξ_y, η_x, η_y) appear **everywhere** in the transformed equations!

---

## **Approach 2: Strong Conservation Form** (Alternative)

Instead of transforming derivatives, transform the entire conservation laws:

```cpp
∂(JQ)/∂t + ∂(JE)/∂ξ + ∂(JF)/∂η = JSource

where:
  J = Jacobian of transformation
  Q = [u, v, p] (conserved variables)
  E, F = flux vectors (include metric terms)
```

---

## 📊 Side-by-Side Comparison

| Aspect | Current Code (Cartesian) | Proper Curvilinear |
|--------|-------------------------|-------------------|
| **Grid in solver** | Pretends it's rectangular | Actually IS rectangular (computational space) |
| **Spacing** | Uses constant `dx = 0.0083 m` | Uses uniform `Δξ = 1`, `Δη = 1` |
| **Derivatives** | `du/dx = (u[j+1]-u[j-1])/(2*dx)` | `du/dx = ξ_x·∂u/∂ξ + η_x·∂u/∂η` |
| **Jacobian terms** | **NONE** ❌ | **Everywhere** ✓ |
| **Metric storage** | No metrics needed | Need ξ_x, ξ_y, η_x, η_y, J at each cell |
| **Accuracy** | Wrong unless grid is Cartesian | Exact for curvilinear grids |

---

## 🔍 What's Actually Happening in Your Code

### Current Workflow (INCORRECT for O-grid):

```
1. Generate O-grid (curved)
   └─> mesh_generator.cpp creates grid.x[i][j], grid.y[i][j]
   
2. Compute "average" spacing
   └─> dx = (xmax - xmin) / nx  ← Single number!
   └─> dy = (ymax - ymin) / ny  ← Doesn't vary per cell!
   
3. Solve Navier-Stokes
   └─> Use dx, dy everywhere (constant)
   └─> Pretend grid is Cartesian
   └─> ❌ Ignores actual grid curvature!
   
4. Results
   └─> Probably "okay" if grid isn't too distorted
   └─> But NOT actually using the O-grid properly!
```

### Proper Curvilinear Workflow (What SHOULD happen):

```
1. Generate O-grid (curved)
   └─> mesh_generator.cpp creates grid.x[i][j], grid.y[i][j]
   
2. Compute metrics everywhere
   └─> For each cell [i][j]:
       ├─> x_ξ = (x[i][j+1] - x[i][j-1]) / (2·Δξ)
       ├─> x_η = (x[i+1][j] - x[i-1][j]) / (2·Δη)
       ├─> y_ξ = (y[i][j+1] - y[i][j-1]) / (2·Δξ)
       ├─> y_η = (y[i+1][j] - y[i-1][j]) / (2·Δη)
       ├─> J = x_ξ·y_η - x_η·y_ξ  (Jacobian)
       ├─> ξ_x = y_η / J   (inverse metrics)
       ├─> ξ_y = -x_η / J
       ├─> η_x = -y_ξ / J
       └─> η_y = x_ξ / J
   
3. Solve Navier-Stokes in computational space (ξ,η)
   └─> Grid is now uniform rectangle!
   └─> Δξ = 1, Δη = 1 everywhere
   └─> But use metrics for all derivatives:
       ∂u/∂x = ξ_x·∂u/∂ξ + η_x·∂u/∂η
       └─────┬─────┘   └─────┬─────┘
          computed       computed
          from metrics   easily!
   
4. Results
   └─> ✓ Properly accounts for grid curvature
   └─> ✓ Exact transformation
   └─> ✓ Actually using the O-grid!
```

---

## 💡 Key Insight: Your Understanding is CORRECT!

You said:
> "do we do the calculations in computational space with like everything being done around the rectangle, pressure being solved around it etc, so yea then the values get transformed with the jacobean multiplier?"

**YES!** That's exactly how it SHOULD work:

1. **Solve in computational space** (perfect rectangle, uniform Δξ and Δη)
2. **Transform derivatives** using Jacobian/metrics
3. **Jacobian appears everywhere** in the equations

---

## 🛠️ What Needs to Change in Your Code

### Current:
```cpp
// Assumes Cartesian grid
double du_dx = (u_right - u_left) / (2.0 * dx);  // ❌
```

### Proper Curvilinear:
```cpp
// Computational derivatives (uniform spacing Δξ = 1)
double du_dxi = (u[i][j+1] - u[i][j-1]) / 2.0;    // ∂u/∂ξ
double du_deta = (u[i+1][j] - u[i-1][j]) / 2.0;   // ∂u/∂η

// Transform to physical using metrics
double du_dx = metrics.xi_x[i][j] * du_dxi + metrics.eta_x[i][j] * du_deta;  // ✓
//             └─────────── Jacobian multipliers! ───────────┘
```

---

## 📐 Example with Numbers

### Current code at cell [30][50]:
```cpp
dx = 0.0083 m  (constant everywhere)
du_dx = (5.5 - 4.9) / (2.0 * 0.0083) = 36.1 [1/s]
```

### Proper curvilinear at same cell [30][50]:
```cpp
// Computational space (uniform)
Δξ = 1.0, Δη = 1.0
du_dxi = (u[30][51] - u[30][49]) / 2.0 = (5.5 - 4.9) / 2.0 = 0.3
du_deta = (u[31][50] - u[29][50]) / 2.0 = (5.3 - 5.1) / 2.0 = 0.1

// Metrics at this cell (computed from grid geometry)
ξ_x[30][50] = 95.2  [1/m]
η_x[30][50] = 12.3  [1/m]

// Transform (chain rule)
du_dx = ξ_x * du_dxi + η_x * du_deta
      = 95.2 * 0.3 + 12.3 * 0.1
      = 28.56 + 1.23
      = 29.8 [1/s]  ← Different from Cartesian approximation!
```

The difference is the **effect of grid curvature**!

---

## 🎯 Summary

| Question | Answer |
|----------|--------|
| **Do we solve in computational space (rectangle)?** | We SHOULD, but currently DON'T |
| **Do we transform with Jacobian multipliers?** | We SHOULD, but currently DON'T |
| **What does current code do?** | Treats O-grid as if it's Cartesian (wrong!) |
| **Is your understanding correct?** | **YES!** You described it perfectly! |
| **What's missing?** | Metric terms (ξ_x, ξ_y, η_x, η_y, J) |

---

## 🔧 To Make This Proper:

1. **Add metrics to Grid class**:
   ```cpp
   double xi_x[ny][nx], xi_y[ny][nx];
   double eta_x[ny][nx], eta_y[ny][nx];
   double J[ny][nx];
   ```

2. **Compute metrics after mesh generation**:
   ```cpp
   void Grid::computeMetrics() { /* ... */ }
   ```

3. **Rewrite solver to use transformed equations**:
   ```cpp
   double du_dx = xi_x * du_dxi + eta_x * du_deta;  // Chain rule
   ```

4. **Update Poisson solver** with metric terms in Laplacian

Your current code is a **Cartesian approximation** applied to a curvilinear grid. It might give reasonable results if the grid isn't too distorted, but it's not the mathematically correct curvilinear formulation you described! 🚀
