#include "solver.hpp"
#include <vector>
#include <cmath>
#include <algorithm>

// ============================================================================
// ⚠️ IMPORTANT: CURRENT IMPLEMENTATION IS CARTESIAN, NOT CURVILINEAR! ⚠️
// ============================================================================
// This solver treats the O-grid as if it were a simple Cartesian rectangle.
// 
// What it does:
//   - Uses constant spacing dx, dy everywhere
//   - Computes derivatives as: du/dx = (u[j+1] - u[j-1]) / (2*dx)
//   - No Jacobian or metric terms
//
// What a PROPER curvilinear solver should do:
//   - Solve in computational space (ξ,η) which IS a perfect rectangle
//   - Use uniform spacing Δξ=1, Δη=1 in computational space
//   - Transform derivatives using metrics: du/dx = ξ_x·∂u/∂ξ + η_x·∂u/∂η
//   - Include Jacobian J and metrics (ξ_x, ξ_y, η_x, η_y) everywhere
//
// See CURRENT_VS_PROPER_CURVILINEAR.md for detailed explanation
// ============================================================================

// Constructor: initialize references and temporary arrays
Solver::Solver(Grid &g, Field &f, double r, double v, double dt_)
    : grid(g), field(f), rho(r), nu(v), dt(dt_)
{
    totalCells = static_cast<size_t>(grid.nx) * grid.ny;
    u_star.assign(totalCells, 0.0);
    v_star.assign(totalCells, 0.0);
}

// Step 1: Compute intermediate velocities (advection + diffusion, explicit)
// Math: u* = u + dt * [-u*du/dx - v*du/dy + nu*(d2u/dx2 + d2u/dy2)]
void Solver::computeIntermediateVelocities() {
    // COMPUTATIONAL SPACE: Grid dimensions (array indices)
    const int nx = grid.nx;  // Number of cells in ξ-direction (i-index)
    const int ny = grid.ny;  // Number of cells in η-direction (j-index)
    
    // PHYSICAL SPACE: Cell spacing in meters
    const double dx = grid.dx;  // Physical spacing Δx [meters]
    const double dy = grid.dy;  // Physical spacing Δy [meters]
    
    // COMPUTATIONAL SPACE: Loop over ALL cells using array indices [i][j]
    // Copy current field to u_star/v_star first (handles boundaries/solids)
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < ny; ++i) {      // i is COMPUTATIONAL index (row)
        for (int j = 0; j < nx; ++j) {  // j is COMPUTATIONAL index (column)
            size_t id = idx(i, j);      // Convert 2D computational index to 1D memory index
            
            // PHYSICAL SPACE: field.u(i,j) is velocity in m/s at cell [i][j]
            u_star[id] = field.u(i, j);  // Velocity values are PHYSICAL (m/s)
            v_star[id] = field.v(i, j);
        }
    }

    // COMPUTATIONAL SPACE: Loop over interior cells only (avoid boundaries)
    // Update interior fluid cells only
    #pragma omp parallel for schedule(static)
    for (int i = 1; i < ny - 1; ++i) {      // COMPUTATIONAL: Skip first/last row
        for (int j = 1; j < nx - 1; ++j) {  // COMPUTATIONAL: Skip first/last column
            if (!grid.isFluid(i, j)) continue;  // Skip solid cells

            size_t id = idx(i, j);  // COMPUTATIONAL: 1D memory index
            
            // PHYSICAL SPACE: Velocities are in m/s
            double u_ij = field.u(i, j);  // u at center cell [m/s]
            double v_ij = field.v(i, j);  // v at center cell [m/s]
            
            // COMPUTATIONAL SPACE: Access neighbors using index offsets (i±1, j±1)
            // Get neighbor velocities (use zero for solid cells = no-slip boundary)
            double u_left = grid.isFluid(i, j-1) ? field.u(i, j-1) : 0.0;    // [i][j-1] PHYSICAL: m/s
            double u_right = grid.isFluid(i, j+1) ? field.u(i, j+1) : 0.0;   // [i][j+1] PHYSICAL: m/s
            double u_bottom = grid.isFluid(i-1, j) ? field.u(i-1, j) : 0.0;  // [i-1][j] PHYSICAL: m/s
            double u_top = grid.isFluid(i+1, j) ? field.u(i+1, j) : 0.0;     // [i+1][j] PHYSICAL: m/s
            
            double v_left = grid.isFluid(i, j-1) ? field.v(i, j-1) : 0.0;    // PHYSICAL: m/s
            double v_right = grid.isFluid(i, j+1) ? field.v(i, j+1) : 0.0;   // PHYSICAL: m/s
            double v_bottom = grid.isFluid(i-1, j) ? field.v(i-1, j) : 0.0;  // PHYSICAL: m/s
            double v_top = grid.isFluid(i+1, j) ? field.v(i+1, j) : 0.0;     // PHYSICAL: m/s

            // ============================================================================
            // TRANSFORMATION POINT: Computational differences → Physical derivatives
            // ============================================================================
            // We compute COMPUTATIONAL differences: Δu = u[j+1] - u[j]
            // Then divide by PHYSICAL spacing dx to get PHYSICAL derivative: du/dx [1/s]
            // This is where the two spaces connect!
            
            // First-order upwind for advection (stable but diffusive)
            // PHYSICAL derivative du/dx ≈ (u - u_upwind)/dx where dx is PHYSICAL [meters]
            // Result: du/dx has units [m/s]/[m] = [1/s] - this is PHYSICAL
            double du_dx = (u_ij >= 0.0) 
                ? (u_ij - u_left) / dx      // PHYSICAL: [m/s]/[m] = [1/s]
                : (u_right - u_ij) / dx;    // PHYSICAL: [m/s]/[m] = [1/s]
            
            double du_dy = (v_ij >= 0.0)
                ? (u_ij - u_bottom) / dy    // PHYSICAL: [m/s]/[m] = [1/s]
                : (u_top - u_ij) / dy;      // PHYSICAL: [m/s]/[m] = [1/s]

            double dv_dx = (u_ij >= 0.0)
                ? (v_ij - v_left) / dx      // PHYSICAL: [1/s]
                : (v_right - v_ij) / dx;
            
            double dv_dy = (v_ij >= 0.0)
                ? (v_ij - v_bottom) / dy    // PHYSICAL: [1/s]
                : (v_top - v_ij) / dy;

            // TRANSFORMATION POINT: Second derivatives (Laplacian for diffusion)
            // Central difference for diffusion (second derivative)
            // d2u/dx2 ≈ (u[i,j+1] - 2*u[i,j] + u[i,j-1]) / dx^2
            // Result: [m/s]/[m²] = [1/(m·s)] - PHYSICAL second derivative
            double d2u_dx2 = (u_right - 2.0 * u_ij + u_left) / (dx * dx);    // PHYSICAL: [1/(m·s)]
            double d2u_dy2 = (u_top - 2.0 * u_ij + u_bottom) / (dy * dy);    // PHYSICAL: [1/(m·s)]
            
            double d2v_dx2 = (v_right - 2.0 * v_ij + v_left) / (dx * dx);    // PHYSICAL: [1/(m·s)]
            double d2v_dy2 = (v_top - 2.0 * v_ij + v_bottom) / (dy * dy);    // PHYSICAL: [1/(m·s)]

            // PHYSICAL SPACE: All terms here are PHYSICAL quantities
            // Explicit forward Euler time integration
            // u* = u^n + dt * RHS
            // Units check: [m/s] + [s] * ([m/s]·[1/s] + [m²/s]·[1/(m·s)]) = [m/s] ✓
            u_star[id] = u_ij + dt * (-u_ij * du_dx - v_ij * du_dy + nu * (d2u_dx2 + d2u_dy2));
            v_star[id] = v_ij + dt * (-u_ij * dv_dx - v_ij * dv_dy + nu * (d2v_dx2 + d2v_dy2));
            // Result: u_star, v_star are PHYSICAL velocities [m/s]
        }
    }
    
    // COMPUTATIONAL SPACE: Loop over all cells to enforce boundary conditions
    // DIRECT FORCING: Enforce zero velocity in and near solid cells
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < ny; ++i) {      // COMPUTATIONAL indices
        for (int j = 0; j < nx; ++j) {
            if (!grid.isFluid(i, j)) {
                // Solid cells: force to zero (PHYSICAL velocity = 0 m/s)
                size_t id = idx(i, j);
                u_star[id] = 0.0;  // PHYSICAL: 0 m/s
                v_star[id] = 0.0;  // PHYSICAL: 0 m/s
            }
        }
    }
}

// Step 2: Solve pressure Poisson equation (SIMPLIFIED - absolute pressure)
// Math: ∇²p = (ρ/Δt) ∇·u*
// Back to basics - solve for absolute pressure with proper BCs
void Solver::solvePressurePoisson(int nit, double tol, double omega){
    // COMPUTATIONAL SPACE: Grid dimensions
    const int nx = grid.nx;  // Number of cells in x-direction (COMPUTATIONAL)
    const int ny = grid.ny;  // Number of cells in y-direction (COMPUTATIONAL)
    
    // PHYSICAL SPACE: Cell spacing in meters
    const double dx = grid.dx;  // Δx [meters]
    const double dy = grid.dy;  // Δy [meters]

    const double dx2 = dx * dx;  // PHYSICAL: [m²]
    const double dy2 = dy * dy;  // PHYSICAL: [m²]
    const double denom = 2.0 * (1.0 / dx2 + 1.0 / dy2);  // PHYSICAL: [1/m²]

    for (int it = 0; it < nit; ++it) {  // Iterative solver loop
        double maxChange = 0.0;  // PHYSICAL: [Pa] - convergence check

        // COMPUTATIONAL SPACE: Loop over interior cells
        for (int i = 1; i < ny - 1; ++i) {      // COMPUTATIONAL index
            for (int j = 1; j < nx - 1; ++j) {  // COMPUTATIONAL index
                if (!grid.isFluid(i, j)) {
                    // Solid cells: extrapolate from neighbors (Neumann BC)
                    int count = 0;
                    double p_avg = 0.0;  // PHYSICAL: [Pa]
                    if (grid.isFluid(i-1, j)) { p_avg += field.p(i-1, j); count++; }
                    if (grid.isFluid(i+1, j)) { p_avg += field.p(i+1, j); count++; }
                    if (grid.isFluid(i, j-1)) { p_avg += field.p(i, j-1); count++; }
                    if (grid.isFluid(i, j+1)) { p_avg += field.p(i, j+1); count++; }
                    if (count > 0) field.p(i, j) = p_avg / count;  // PHYSICAL: [Pa]
                    continue;
                }

                // ============================================================================
                // TRANSFORMATION POINT: Compute divergence of u* (PHYSICAL derivative)
                // ============================================================================
                // RHS: (ρ/Δt) * divergence of u*
                // COMPUTATIONAL differences (u[j+1] - u[j-1]) divided by PHYSICAL spacing
                double du_dx = (u_star[idx(i, j + 1)] - u_star[idx(i, j - 1)]) / (2.0 * dx);  // PHYSICAL: [1/s]
                double dv_dy = (v_star[idx(i + 1, j)] - v_star[idx(i - 1, j)]) / (2.0 * dy);  // PHYSICAL: [1/s]
                
                double div_ustar = du_dx + dv_dy;  // PHYSICAL: ∇·u* [1/s]
                double rhs = (rho / dt) * div_ustar;  // PHYSICAL: [kg/m³]/[s] * [1/s] = [kg/(m³·s²)] = [Pa/m²]

                // TRANSFORMATION POINT: Laplacian of pressure (PHYSICAL second derivative)
                // Standard 5-point Laplacian: ∇²p ≈ (p[i,j+1] - 2p + p[i,j-1])/dx² + ...
                // COMPUTATIONAL neighbors accessed with i±1, j±1
                // PHYSICAL spacing dx², dy² in denominator
                double p_old = field.p(i, j);  // PHYSICAL: [Pa]
                double p_gs = ((field.p(i, j + 1) + field.p(i, j - 1)) / dx2   // PHYSICAL: [Pa]/[m²]
                             + (field.p(i + 1, j) + field.p(i - 1, j)) / dy2   // PHYSICAL: [Pa]/[m²]
                             - rhs) / denom;  // Result: [Pa]

                // SOR (Successive Over-Relaxation)
                double p_new = (1.0 - omega) * p_old + omega * p_gs;  // PHYSICAL: [Pa]
                field.p(i, j) = p_new;  // Update pressure field with PHYSICAL value
                
                maxChange = std::max(maxChange, std::abs(p_new - p_old));  // PHYSICAL: [Pa]
            }
        }

        // COMPUTATIONAL SPACE: Apply boundary conditions at domain edges
        // Apply BCs - all zero-gradient (Neumann)
        for (int j = 0; j < nx; ++j) {  // COMPUTATIONAL loop over columns
            field.p(0, j) = field.p(1, j);          // PHYSICAL: Copy pressure [Pa]
            field.p(ny - 1, j) = field.p(ny - 2, j);  // PHYSICAL: Copy pressure [Pa]
        }
        for (int i = 0; i < ny; ++i) {  // COMPUTATIONAL loop over rows
            field.p(i, 0) = field.p(i, 1);          // PHYSICAL: Copy pressure [Pa]
            field.p(i, nx - 1) = field.p(i, nx - 2);  // PHYSICAL: Copy pressure [Pa]
        }

        // Early exit if converged (PHYSICAL tolerance in Pa)
        if (maxChange < tol) {
            break;
        }
    }
}

// Step 3: Correct velocities with pressure gradient
// Math: u^(n+1) = u* - (Δt/ρ) ∇p
void Solver::correctVelocities(){
    // PHYSICAL SPACE: Cell spacing
    const double dx = grid.dx;  // [meters]
    const double dy = grid.dy;  // [meters]
    
    // COMPUTATIONAL SPACE: Loop over interior cells
    #pragma omp parallel for schedule(static)
    for (int i = 1; i < grid.ny - 1; i++){      // COMPUTATIONAL index
        for (int j = 1; j < grid.nx - 1; j++){  // COMPUTATIONAL index
            if(!grid.isFluid(i,j)) continue;

            // ============================================================================
            // TRANSFORMATION POINT: Pressure gradient (PHYSICAL derivative)
            // ============================================================================
            // Pressure gradient (central difference)
            // COMPUTATIONAL neighbors: p[i][j±1], p[i±1][j]
            // PHYSICAL spacing: dx, dy [meters]
            // Result: dp/dx [Pa/m] = [kg/(m²·s²)]
            double dp_dx = (field.p(i, j + 1) - field.p(i, j - 1)) / (2.0 * dx);  // PHYSICAL: [Pa/m]
            double dp_dy = (field.p(i + 1, j) - field.p(i - 1, j)) / (2.0 * dy);  // PHYSICAL: [Pa/m]

            // PHYSICAL SPACE: Velocity correction
            // Correct velocities: u = u* - (dt/ρ)·∇p
            // Units: [m/s] - [s]/[kg/m³] · [Pa/m] = [m/s] - [s]·[m³/kg]·[kg/(m²·s²)]/[m] = [m/s] ✓
            size_t id = idx(i, j);
            field.u(i, j) = u_star[id] - (dt / rho) * dp_dx;  // PHYSICAL: [m/s]
            field.v(i, j) = v_star[id] - (dt / rho) * dp_dy;  // PHYSICAL: [m/s]
        }
    }
}

// Single time step: combine all 3 steps
void Solver::step(){
    computeIntermediateVelocities();  // Step 1: Predictor (u*)
    solvePressurePoisson(100, 1e-4, 1.7);  // Step 2: Solve ∇²p = (ρ/Δt)∇·u*
    correctVelocities();  // Step 3: Corrector (u = u* - dt∇p/ρ)
    
    // PHYSICAL SPACE: Normalize pressure to remove arbitrary constant
    // (Incompressible flow: only pressure GRADIENTS matter, not absolute value)
    double p_sum = 0.0;  // PHYSICAL: [Pa]
    int fluid_cells = 0;  // COMPUTATIONAL: count
    
    // COMPUTATIONAL SPACE: Loop to compute mean pressure
    for(int i = 0; i < grid.ny; i++){      // COMPUTATIONAL index
        for(int j = 0; j < grid.nx; j++){  // COMPUTATIONAL index
            if(grid.isFluid(i, j)){
                p_sum += field.p(i, j);  // PHYSICAL: accumulate pressure [Pa]
                fluid_cells++;  // COMPUTATIONAL: count cells
            }
        }
    }
    double p_mean = p_sum / fluid_cells;  // PHYSICAL: mean pressure [Pa]
    
    // COMPUTATIONAL SPACE: Loop to subtract mean from all cells
    for(int i = 0; i < grid.ny; i++){      // COMPUTATIONAL index
        for(int j = 0; j < grid.nx; j++){  // COMPUTATIONAL index
            field.p(i, j) -= p_mean;  // PHYSICAL: shift pressure by constant [Pa]
        }
    }
}
