# CFD Solver - Phase 1 Summary & Phase 2 Planning

**Date**: November 9-10, 2025  
**Status**: Phase 1 Complete ✅ | Phase 2 Planning 📋

---

## Phase 1: SUCCESS - Cylinder Vortex Shedding

### What Worked ✅

1. **Fractional-Step Method Implementation**
   - 3-step projection method (advection-diffusion → Poisson → correction)
   - Explicit forward Euler time integration
   - SOR Poisson solver (100 iterations, ω=1.7)

2. **Circular Cylinder Test Case**
   - **Reynolds number**: Re = 1000 (U∞=5 m/s, D=0.2m, ν=1e-3)
   - **Grid**: 400×200 cells, dx=dy=0.01m
   - **Results**: Beautiful von Kármán vortex street!
   - **Physics**: Clear vortex shedding, Strouhal number ~0.2
   - **Stability**: 20,000 timesteps with no NaNs, mass conserved

3. **Visualization**
   - GIF animation of vortex evolution
   - Velocity statistics showing oscillations
   - Clean flow patterns

### Key Achievements
- ✅ Stable incompressible flow solver
- ✅ Proper vortex shedding dynamics
- ✅ OpenMP parallelization working
- ✅ NaN-free solid boundary treatment
- ✅ Mass conservation maintained

---

## Phase 1: PARTIAL SUCCESS - Airfoil Flow

### What We Tried

1. **NACA 2412 Airfoil**
   - Chord = 0.6m, α = -5°, Re = 15,000
   - Grid: 400×200 cells
   - Solid cells: 303 (properly masked)

2. **Boundary Condition Attempts**
   - ❌ Zero-gradient BCs everywhere → singular Poisson system
   - ❌ Outlet p=0 Dirichlet → artificial pressure gradient
   - ❌ Pressure normalization by mean → still has left-right gradient
   - ✅ Direct forcing (u*=0, v*=0 in solids) → stable but physics wrong

### Problems Identified

1. **Pressure Field Issues**
   - Large artificial gradient from inlet (-8 Pa) to outlet (+4 Pa)
   - Pressure coefficient doesn't show airfoil suction/pressure surfaces
   - Cp values dominated by domain-wide gradient, not local flow

2. **Velocity Field Issues**
   - Flow decelerates around airfoil instead of accelerating
   - No clear acceleration over top surface
   - Wake behavior incorrect
   - Centerline velocity shows bizarre jumps

3. **Root Cause Analysis**
   - **Poisson solver BCs are wrong**: Zero-gradient on all boundaries creates mathematically singular system
   - **Pressure reference issue**: Normalizing by mean doesn't fix the gradient structure
   - **Immersed boundary coupling**: Solid cells still affecting fluid pressure field incorrectly

---

## Phase 2: Boundary Condition Fix Strategy

### The Fundamental Problem

**Current situation:**
```
Boundaries: All Neumann (∂p/∂n = 0)
Result: ∇²p = f has infinite solutions (p + constant)
Our fix: Subtract mean pressure → removes constant but gradient structure remains wrong
```

**What we need:**
- One Dirichlet BC (p = value) to anchor the pressure
- But it must be done WITHOUT creating artificial gradients

### Proposed Solutions (Ranked by Difficulty)

---

#### **Option 1: Proper Outlet Pressure BC** ⭐ RECOMMENDED
**Difficulty**: Easy  
**Implementation time**: 30 min  

**Strategy:**
1. Set **outlet pressure to zero** as Dirichlet BC: `p(i, nx-1) = 0.0`
2. BUT use **extrapolation** for interior cells near outlet instead of Poisson solve
3. Apply Poisson only for `j = 1 to nx-2`, skip `j = nx-1`
4. This gives proper reference WITHOUT letting the BC contaminate the interior

**Math:**
```cpp
// In Poisson solver:
for (j = 1; j < nx-1; ++j) {  // Stop before outlet!
    // Solve Poisson normally
}
// Then set outlet explicitly:
for (i = 0; i < ny; ++i) {
    field.p(i, nx-1) = 0.0;  // Dirichlet
}
```

**Expected result:**
- Pressure anchored to zero at outlet
- Interior pressure determined by flow physics
- No artificial gradient through domain

---

#### **Option 2: Incremental Pressure Correction**
**Difficulty**: Medium  
**Implementation time**: 2 hours  

**Strategy:**
Solve for **pressure correction** φ instead of absolute pressure p:
```
Step 1: u* = u^n + dt*[advection + diffusion]
Step 2: ∇²φ = (ρ/dt)∇·u*
Step 3: u^(n+1) = u* - (dt/ρ)∇φ
        p^(n+1) = p^n + φ
```

**Advantages:**
- φ has well-defined BCs: φ=0 on all boundaries
- No pressure drift issues
- Standard in modern CFD codes

**Changes needed:**
- Add pressure correction field `phi`
- Modify step() function to use incremental approach
- Store previous pressure

---

#### **Option 3: Convective Outlet BC**
**Difficulty**: Medium-Hard  
**Implementation time**: 3-4 hours  

**Strategy:**
Apply **convective** (wave) BC at outlet:
```
∂p/∂t + U_conv * ∂p/∂x = 0
```
Discretized as:
```cpp
p(i,nx-1) = p_old(i,nx-1) - U_conv*dt/dx * (p(i,nx-1) - p(i,nx-2))
```

**Advantages:**
- Allows disturbances to leave domain smoothly
- Good for unsteady flows
- No artificial reflection

**Disadvantages:**
- Needs to store `p_old`
- Choice of `U_conv` is empirical (usually U_infty)

---

#### **Option 4: Body-Fitted Grid** 
**Difficulty**: HARD 🔥  
**Implementation time**: Weeks  

**Why we're NOT doing this:**
- Need mesh generator
- Curved coordinate transformations
- Different mesh for each geometry
- Way too complex for learning project

---

#### **Option 5: Advanced IB Method (Feedback Forcing)**
**Difficulty**: HARD 🔥  
**Implementation time**: Days  

**What it is:**
Add forcing term to momentum equation:
```
∂u/∂t + u·∇u = -∇p/ρ + ν∇²u + f_IB
```
where `f_IB` is computed to enforce u=0 inside solid

**Why later:**
- Requires iterative force calculation
- Needs careful tuning
- More of a research topic

---

## Recommended Path Forward

### **Phase 2.1: Fix Outlet Pressure BC** (Next Session)

**Task list:**
1. ✅ Document current status (this file!)
2. ⬜ Implement Option 1 (proper outlet Dirichlet BC)
3. ⬜ Test with cylinder (should still work)
4. ⬜ Test with airfoil (should see proper pressure field!)
5. ⬜ If successful: Calculate lift/drag coefficients
6. ⬜ Compare with theory (NACA 2412 at α=5° should give Cl ~ 0.6-0.8)

**Success criteria:**
- Pressure field shows suction on top, pressure on bottom
- Cp < 0 on top surface, Cp > 0 on bottom surface  
- No artificial left-right gradient
- Flow acceleration over top, deceleration below

---

### **Phase 2.2: Validation** (After BC fix)

Once pressure is working:

1. **Airfoil validation**
   - Run NACA 0012 at α=0° (should have Cl=0 by symmetry)
   - Run NACA 2412 at α=5° (should have Cl~0.7)
   - Check if Cp distribution matches theory

2. **Grid refinement study**
   - Try 800×400 grid
   - Check if solution converges

3. **Reynolds number sweep**
   - Test Re = 5k, 10k, 20k
   - Watch for flow separation

---

## Technical Lessons Learned

### What Makes CFD Hard

1. **Boundary Conditions Are Everything**
   - Wrong BCs → wrong physics, even if solver is stable
   - Poisson equation needs at least one Dirichlet BC
   - BC choice affects entire solution structure

2. **Immersed Boundary Methods Are Tricky**
   - Simple masking works for thick, simple shapes (cylinder ✅)
   - Thin, complex shapes need more sophisticated treatment
   - Solid cells affect pressure field in subtle ways

3. **Stability ≠ Accuracy**
   - Our solver is stable (no NaNs, mass conserved)
   - But pressure field is still wrong (artificial gradients)
   - Need both stability AND correct physics!

---

## Code Status

### What's Currently Implemented

**Solver (`src/solver.cpp`):**
- ✅ Fractional-step projection method
- ✅ Explicit advection-diffusion with upwind
- ✅ SOR Poisson solver
- ✅ Velocity correction step
- ✅ Direct forcing for solid cells
- ✅ Pressure normalization (mean = 0)
- ❌ Proper outlet pressure BC (NEXT!)

**Geometry (`src/main.cpp` + `src/grid.cpp`):**
- ✅ NACA 4-digit airfoil generator
- ✅ Rotation support
- ✅ Solid cell masking

**Boundary Conditions:**
- ✅ Inlet: Dirichlet (u=U∞, v=0)
- ✅ Top/Bottom: Free-slip
- ✅ Outlet: Zero-gradient (u,v)
- ❌ Outlet pressure: Currently zero-gradient (WRONG!)

**Output:**
- ✅ u, v, p fields saved as CSV
- ✅ NaN masking for visualization
- ✅ Python analysis scripts

---

## Next Session Plan

### Immediate Goal: Fix Outlet Pressure BC

**File to modify**: `src/solver.cpp` - `solvePressurePoisson()`

**Changes needed**:
```cpp
// OLD (current):
for (i = 0; i < ny; ++i) {
    field.p(i, nx-1) = field.p(i, nx-2);  // Zero-gradient
}

// NEW (Option 1):
// Don't solve Poisson at outlet, just set it:
for (int j = 1; j < nx-1; ++j) {  // STOP BEFORE OUTLET
    // ... existing Poisson solver code ...
}
// After solver loop, set outlet:
for (int i = 0; i < ny; ++i) {
    field.p(i, nx-1) = 0.0;  // Dirichlet reference
}
```

**Test procedure**:
1. Run cylinder case → should still work
2. Run airfoil case → check pressure field
3. If pressure gradient disappears → SUCCESS!
4. Calculate Cp and check for airfoil suction peak

---

## Long-term Goals (Phase 3+)

1. **GPU Acceleration** (user mentioned this!)
   - CUDA kernel for Poisson solver
   - GPU-friendly data structures
   - Could get 10-100× speedup

2. **Higher-Order Methods**
   - RK4 time integration (already tried, removed for speed)
   - Higher-order spatial schemes (QUICK, MUSCL)
   - Better accuracy without finer grids

3. **Turbulence Modeling**
   - LES (Large Eddy Simulation)
   - Allow higher Reynolds numbers
   - See turbulent vortex shedding!

4. **More Complex Geometries**
   - Multi-element airfoils
   - 3D wings
   - Moving bodies (flapping, rotating)

---

## Success Metrics

### Phase 1 (Complete) ✅
- [x] Stable solver (no NaNs)
- [x] Mass conservation
- [x] Cylinder vortex shedding
- [x] Visualization working
- [x] OpenMP parallelization

### Phase 2 (In Progress) 🔄
- [ ] Proper pressure field around airfoil
- [ ] Cp distribution matches theory
- [ ] Lift coefficient calculation
- [ ] Flow acceleration on suction side

### Phase 3 (Future) 📋
- [ ] GPU acceleration
- [ ] Turbulence modeling
- [ ] Higher Re simulations
- [ ] Complex geometries

---

## References for Phase 2

### Textbooks
1. **Ferziger & Perić** - "Computational Methods for Fluid Dynamics"
   - Chapter on pressure BCs (Section 7.3)
   - Fractional-step methods (Section 6.3)

2. **Tannehill et al.** - "Computational Fluid Mechanics and Heat Transfer"
   - Immersed boundary methods

### Papers
1. **Mittal & Iaccarino (2005)** - "Immersed boundary methods" (Annual Review)
   - Direct forcing method (what we implemented)
   - Feedback forcing (more advanced)

2. **Kim & Moin (1985)** - "Application of fractional-step method to incompressible flows"
   - Classic projection method paper
   - Pressure BC treatment

---

**Status**: Phase 1 complete, ready for Phase 2.1! 🚀
