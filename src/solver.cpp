#include "solver.hpp"
#include <vector>
#include <cmath>
#include <algorithm>

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
    const int nx = grid.nx;
    const int ny = grid.ny;
    const double dx = grid.dx;
    const double dy = grid.dy;
    
    // Copy current field to u_star/v_star first (handles boundaries/solids)
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < ny; ++i) {
        for (int j = 0; j < nx; ++j) {
            size_t id = idx(i, j);
            u_star[id] = field.u(i, j);
            v_star[id] = field.v(i, j);
        }
    }

    // Update interior fluid cells only
    #pragma omp parallel for schedule(static)
    for (int i = 1; i < ny - 1; ++i) {
        for (int j = 1; j < nx - 1; ++j) {
            if (!grid.isFluid(i, j)) continue;

            size_t id = idx(i, j);
            double u_ij = field.u(i, j);
            double v_ij = field.v(i, j);
            
            // Get neighbor velocities (use zero for solid cells = no-slip boundary)
            double u_left = grid.isFluid(i, j-1) ? field.u(i, j-1) : 0.0;
            double u_right = grid.isFluid(i, j+1) ? field.u(i, j+1) : 0.0;
            double u_bottom = grid.isFluid(i-1, j) ? field.u(i-1, j) : 0.0;
            double u_top = grid.isFluid(i+1, j) ? field.u(i+1, j) : 0.0;
            
            double v_left = grid.isFluid(i, j-1) ? field.v(i, j-1) : 0.0;
            double v_right = grid.isFluid(i, j+1) ? field.v(i, j+1) : 0.0;
            double v_bottom = grid.isFluid(i-1, j) ? field.v(i-1, j) : 0.0;
            double v_top = grid.isFluid(i+1, j) ? field.v(i+1, j) : 0.0;

            // First-order upwind for advection (stable but diffusive)
            // du/dx ≈ (u - u_upwind)/dx where upwind is based on velocity direction
            double du_dx = (u_ij >= 0.0) 
                ? (u_ij - u_left) / dx
                : (u_right - u_ij) / dx;
            
            double du_dy = (v_ij >= 0.0)
                ? (u_ij - u_bottom) / dy
                : (u_top - u_ij) / dy;

            double dv_dx = (u_ij >= 0.0)
                ? (v_ij - v_left) / dx
                : (v_right - v_ij) / dx;
            
            double dv_dy = (v_ij >= 0.0)
                ? (v_ij - v_bottom) / dy
                : (v_top - v_ij) / dy;

            // Central difference for diffusion (second derivative)
            // d2u/dx2 ≈ (u[i,j+1] - 2*u[i,j] + u[i,j-1]) / dx^2
            double d2u_dx2 = (u_right - 2.0 * u_ij + u_left) / (dx * dx);
            double d2u_dy2 = (u_top - 2.0 * u_ij + u_bottom) / (dy * dy);
            
            double d2v_dx2 = (v_right - 2.0 * v_ij + v_left) / (dx * dx);
            double d2v_dy2 = (v_top - 2.0 * v_ij + v_bottom) / (dy * dy);

            // Explicit forward Euler time integration
            // u* = u^n + dt * RHS
            u_star[id] = u_ij + dt * (-u_ij * du_dx - v_ij * du_dy + nu * (d2u_dx2 + d2u_dy2));
            v_star[id] = v_ij + dt * (-u_ij * dv_dx - v_ij * dv_dy + nu * (d2v_dx2 + d2v_dy2));
        }
    }
    
    // DIRECT FORCING: Enforce zero velocity in and near solid cells
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < ny; ++i) {
        for (int j = 0; j < nx; ++j) {
            if (!grid.isFluid(i, j)) {
                // Solid cells: force to zero
                size_t id = idx(i, j);
                u_star[id] = 0.0;
                v_star[id] = 0.0;
            }
        }
    }
}

// Step 2: Solve pressure Poisson equation (SIMPLIFIED - absolute pressure)
// Math: ∇²p = (ρ/Δt) ∇·u*
// Back to basics - solve for absolute pressure with proper BCs
void Solver::solvePressurePoisson(int nit, double tol, double omega){
    const int nx = grid.nx;
    const int ny = grid.ny;
    const double dx = grid.dx;
    const double dy = grid.dy;

    const double dx2 = dx * dx;
    const double dy2 = dy * dy;
    const double denom = 2.0 * (1.0 / dx2 + 1.0 / dy2);

    for (int it = 0; it < nit; ++it) {
        double maxChange = 0.0;

        for (int i = 1; i < ny - 1; ++i) {
            for (int j = 1; j < nx - 1; ++j) {
                if (!grid.isFluid(i, j)) {
                    // Solid cells: extrapolate from neighbors
                    int count = 0;
                    double p_avg = 0.0;
                    if (grid.isFluid(i-1, j)) { p_avg += field.p(i-1, j); count++; }
                    if (grid.isFluid(i+1, j)) { p_avg += field.p(i+1, j); count++; }
                    if (grid.isFluid(i, j-1)) { p_avg += field.p(i, j-1); count++; }
                    if (grid.isFluid(i, j+1)) { p_avg += field.p(i, j+1); count++; }
                    if (count > 0) field.p(i, j) = p_avg / count;
                    continue;
                }

                // RHS: (ρ/Δt) * divergence of u*
                double du_dx = (u_star[idx(i, j + 1)] - u_star[idx(i, j - 1)]) / (2.0 * dx);
                double dv_dy = (v_star[idx(i + 1, j)] - v_star[idx(i - 1, j)]) / (2.0 * dy);
                
                double div_ustar = du_dx + dv_dy;
                double rhs = (rho / dt) * div_ustar;

                // Standard 5-point Laplacian
                double p_old = field.p(i, j);
                double p_gs = ((field.p(i, j + 1) + field.p(i, j - 1)) / dx2 
                             + (field.p(i + 1, j) + field.p(i - 1, j)) / dy2 
                             - rhs) / denom;

                // SOR
                double p_new = (1.0 - omega) * p_old + omega * p_gs;
                field.p(i, j) = p_new;
                
                maxChange = std::max(maxChange, std::abs(p_new - p_old));
            }
        }

        // Apply BCs - all zero-gradient
        for (int j = 0; j < nx; ++j) {
            field.p(0, j) = field.p(1, j);
            field.p(ny - 1, j) = field.p(ny - 2, j);
        }
        for (int i = 0; i < ny; ++i) {
            field.p(i, 0) = field.p(i, 1);
            field.p(i, nx - 1) = field.p(i, nx - 2);
        }

        // Early exit if converged
        if (maxChange < tol) {
            break;
        }
    }
}

// Step 3: Correct velocities with pressure gradient
// Math: u^(n+1) = u* - (Δt/ρ) ∇p
void Solver::correctVelocities(){
    const double dx = grid.dx;
    const double dy = grid.dy;
    
    #pragma omp parallel for schedule(static)
    for (int i = 1; i < grid.ny - 1; i++){
        for (int j = 1; j < grid.nx - 1; j++){
            if(!grid.isFluid(i,j)) continue;

            // Pressure gradient (central difference)
            double dp_dx = (field.p(i, j + 1) - field.p(i, j - 1)) / (2.0 * dx);
            double dp_dy = (field.p(i + 1, j) - field.p(i - 1, j)) / (2.0 * dy);

            // Correct velocities
            size_t id = idx(i, j);
            field.u(i, j) = u_star[id] - (dt / rho) * dp_dx;
            field.v(i, j) = v_star[id] - (dt / rho) * dp_dy;
        }
    }
}

// Single time step: combine all 3 steps
void Solver::step(){
    computeIntermediateVelocities();
    solvePressurePoisson(100, 1e-4, 1.7);
    correctVelocities();
    
    // Normalize pressure to remove arbitrary constant
    double p_sum = 0.0;
    int fluid_cells = 0;
    for(int i = 0; i < grid.ny; i++){
        for(int j = 0; j < grid.nx; j++){
            if(grid.isFluid(i, j)){
                p_sum += field.p(i, j);
                fluid_cells++;
            }
        }
    }
    double p_mean = p_sum / fluid_cells;
    
    for(int i = 0; i < grid.ny; i++){
        for(int j = 0; j < grid.nx; j++){
            field.p(i, j) -= p_mean;
        }
    }
}
