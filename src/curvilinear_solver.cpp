#include "curvilinear_solver.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

CurvilinearSolver::CurvilinearSolver(const CurvilinearGrid& grid, double Re, double dt)
    : grid_(grid), Re_(Re), nu_(1.0/Re), dt_(dt) {
    
    // Computational grid is UNIFORM (this is the beauty of body-fitted grids!)
    dxi_ = 1.0 / (grid_.nxi() - 1);
    deta_ = 1.0 / (grid_.neta() - 1);
    
    std::cout << "Curvilinear solver initialized:\n";
    std::cout << "  Grid: " << grid_.nxi() << " x " << grid_.neta() << "\n";
    std::cout << "  Re = " << Re_ << ", nu = " << nu_ << "\n";
    std::cout << "  dt = " << dt_ << "\n";
    std::cout << "  Computational spacing: dξ = " << dxi_ << ", dη = " << deta_ << "\n\n";
}

// CRITICAL FUNCTION: Transform ∂u/∂x from physical to computational space
// 
// THEORY (Thompson Eq. 5.12):
// ∂u/∂x = (∂u/∂ξ)·(∂ξ/∂x) + (∂u/∂η)·(∂η/∂x)
//       = u_ξ · ξ_x + u_η · η_x
//
// WHERE:
// - u_ξ, u_η: computed with central differences in (ξ,η) - EASY!
// - ξ_x, η_x: transformation metrics from grid generator
//
double CurvilinearSolver::du_dx(const CurvilinearField& field, int i, int j) const {
    // Compute ∂u/∂ξ using central difference (standard finite difference!)
    double u_xi = (field.u(i+1, j) - field.u(i-1, j)) / (2.0 * dxi_);
    
    // Compute ∂u/∂η using central difference
    double u_eta = (field.u(i, j+1) - field.u(i, j-1)) / (2.0 * deta_);
    
    // Transform to physical space using metrics (Thompson Eq. 5.12)
    return grid_.xi_x(i, j) * u_xi + grid_.eta_x(i, j) * u_eta;
}

double CurvilinearSolver::du_dy(const CurvilinearField& field, int i, int j) const {
    double u_xi = (field.u(i+1, j) - field.u(i-1, j)) / (2.0 * dxi_);
    double u_eta = (field.u(i, j+1) - field.u(i, j-1)) / (2.0 * deta_);
    
    // ∂u/∂y = u_ξ · ξ_y + u_η · η_y
    return grid_.xi_y(i, j) * u_xi + grid_.eta_y(i, j) * u_eta;
}

double CurvilinearSolver::dv_dx(const CurvilinearField& field, int i, int j) const {
    double v_xi = (field.v(i+1, j) - field.v(i-1, j)) / (2.0 * dxi_);
    double v_eta = (field.v(i, j+1) - field.v(i, j-1)) / (2.0 * deta_);
    return grid_.xi_x(i, j) * v_xi + grid_.eta_x(i, j) * v_eta;
}

double CurvilinearSolver::dv_dy(const CurvilinearField& field, int i, int j) const {
    double v_xi = (field.v(i+1, j) - field.v(i-1, j)) / (2.0 * dxi_);
    double v_eta = (field.v(i, j+1) - field.v(i, j-1)) / (2.0 * deta_);
    return grid_.xi_y(i, j) * v_xi + grid_.eta_y(i, j) * v_eta;
}

// LAPLACIAN TRANSFORMATION
// This is more complex - see Thompson Eq. 5.15-5.16
//
// FULL FORMULA:
// ∇²u = α·u_ξξ + 2β·u_ξη + γ·u_ηη + (metric derivative terms)
//
// where α, β, γ are the SAME as in grid generation:
// α = ξ_x² + ξ_y²
// β = ξ_x·η_x + ξ_y·η_y  
// γ = η_x² + η_y²
//
double CurvilinearSolver::laplacian_u(const CurvilinearField& field, int i, int j) const {
    // Second derivatives in computational space (standard stencils!)
    double u_xixi = (field.u(i+1,j) - 2.0*field.u(i,j) + field.u(i-1,j)) / (dxi_*dxi_);
    double u_etaeta = (field.u(i,j+1) - 2.0*field.u(i,j) + field.u(i,j-1)) / (deta_*deta_);
    double u_xieta = (field.u(i+1,j+1) - field.u(i+1,j-1) 
                    - field.u(i-1,j+1) + field.u(i-1,j-1)) / (4.0*dxi_*deta_);
    
    // Get metrics at this point
    double xi_x = grid_.xi_x(i, j);
    double xi_y = grid_.xi_y(i, j);
    double eta_x = grid_.eta_x(i, j);
    double eta_y = grid_.eta_y(i, j);
    
    // Metric coefficients (Thompson Eq. 5.6)
    double alpha = xi_x*xi_x + xi_y*xi_y;
    double beta = xi_x*eta_x + xi_y*eta_y;
    double gamma = eta_x*eta_x + eta_y*eta_y;
    
    // Clamp extreme metric coefficients to prevent numerical overflow
    // This is needed for highly stretched grids near boundaries
    const double MAX_ALPHA = 1000.0;  // Limit transformation coefficient magnitude
    if (alpha > MAX_ALPHA) alpha = MAX_ALPHA;
    if (gamma > MAX_ALPHA) gamma = MAX_ALPHA;
    if (std::abs(beta) > MAX_ALPHA) beta = (beta > 0) ? MAX_ALPHA : -MAX_ALPHA;
    
    // Transformed Laplacian (Thompson Eq. 5.15)
    // Simplified version (neglecting metric derivative terms for now - small contribution)
    return alpha * u_xixi + 2.0 * beta * u_xieta + gamma * u_etaeta;
}

double CurvilinearSolver::laplacian_v(const CurvilinearField& field, int i, int j) const {
    double v_xixi = (field.v(i+1,j) - 2.0*field.v(i,j) + field.v(i-1,j)) / (dxi_*dxi_);
    double v_etaeta = (field.v(i,j+1) - 2.0*field.v(i,j) + field.v(i,j-1)) / (deta_*deta_);
    double v_xieta = (field.v(i+1,j+1) - field.v(i+1,j-1) 
                    - field.v(i-1,j+1) + field.v(i-1,j-1)) / (4.0*dxi_*deta_);
    
    double xi_x = grid_.xi_x(i, j);
    double xi_y = grid_.xi_y(i, j);
    double eta_x = grid_.eta_x(i, j);
    double eta_y = grid_.eta_y(i, j);
    
    double alpha = xi_x*xi_x + xi_y*xi_y;
    double beta = xi_x*eta_x + xi_y*eta_y;
    double gamma = eta_x*eta_x + eta_y*eta_y;
    
    // Clamp extreme metrics (same as laplacian_u)
    const double MAX_ALPHA = 1000.0;
    if (alpha > MAX_ALPHA) alpha = MAX_ALPHA;
    if (gamma > MAX_ALPHA) gamma = MAX_ALPHA;
    if (std::abs(beta) > MAX_ALPHA) beta = (beta > 0) ? MAX_ALPHA : -MAX_ALPHA;
    
    return alpha * v_xixi + 2.0 * beta * v_xieta + gamma * v_etaeta;
}

// Laplacian of u_star (for implicit diffusion solver)
double CurvilinearSolver::laplacian_u_star(const CurvilinearField& field, int i, int j) const {
    double u_xixi = (field.u_star(i+1,j) - 2.0*field.u_star(i,j) + field.u_star(i-1,j)) / (dxi_*dxi_);
    double u_etaeta = (field.u_star(i,j+1) - 2.0*field.u_star(i,j) + field.u_star(i,j-1)) / (deta_*deta_);
    double u_xieta = (field.u_star(i+1,j+1) - field.u_star(i+1,j-1) 
                    - field.u_star(i-1,j+1) + field.u_star(i-1,j-1)) / (4.0*dxi_*deta_);
    
    double xi_x = grid_.xi_x(i, j);
    double xi_y = grid_.xi_y(i, j);
    double eta_x = grid_.eta_x(i, j);
    double eta_y = grid_.eta_y(i, j);
    
    double alpha = xi_x*xi_x + xi_y*xi_y;
    double beta = xi_x*eta_x + xi_y*eta_y;
    double gamma = eta_x*eta_x + eta_y*eta_y;
    
    const double MAX_ALPHA = 1000.0;
    if (alpha > MAX_ALPHA) alpha = MAX_ALPHA;
    if (gamma > MAX_ALPHA) gamma = MAX_ALPHA;
    if (std::abs(beta) > MAX_ALPHA) beta = (beta > 0) ? MAX_ALPHA : -MAX_ALPHA;
    
    return alpha * u_xixi + 2.0 * beta * u_xieta + gamma * u_etaeta;
}

double CurvilinearSolver::laplacian_v_star(const CurvilinearField& field, int i, int j) const {
    double v_xixi = (field.v_star(i+1,j) - 2.0*field.v_star(i,j) + field.v_star(i-1,j)) / (dxi_*dxi_);
    double v_etaeta = (field.v_star(i,j+1) - 2.0*field.v_star(i,j) + field.v_star(i,j-1)) / (deta_*deta_);
    double v_xieta = (field.v_star(i+1,j+1) - field.v_star(i+1,j-1) 
                    - field.v_star(i-1,j+1) + field.v_star(i-1,j-1)) / (4.0*dxi_*deta_);
    
    double xi_x = grid_.xi_x(i, j);
    double xi_y = grid_.xi_y(i, j);
    double eta_x = grid_.eta_x(i, j);
    double eta_y = grid_.eta_y(i, j);
    
    double alpha = xi_x*xi_x + xi_y*xi_y;
    double beta = xi_x*eta_x + xi_y*eta_y;
    double gamma = eta_x*eta_x + eta_y*eta_y;
    
    const double MAX_ALPHA = 1000.0;
    if (alpha > MAX_ALPHA) alpha = MAX_ALPHA;
    if (gamma > MAX_ALPHA) gamma = MAX_ALPHA;
    if (std::abs(beta) > MAX_ALPHA) beta = (beta > 0) ? MAX_ALPHA : -MAX_ALPHA;
    
    return alpha * v_xixi + 2.0 * beta * v_xieta + gamma * v_etaeta;
}

double CurvilinearSolver::dp_dx(const CurvilinearField& field, int i, int j) const {
    double p_xi = (field.p(i+1, j) - field.p(i-1, j)) / (2.0 * dxi_);
    double p_eta = (field.p(i, j+1) - field.p(i, j-1)) / (2.0 * deta_);
    return grid_.xi_x(i, j) * p_xi + grid_.eta_x(i, j) * p_eta;
}

double CurvilinearSolver::dp_dy(const CurvilinearField& field, int i, int j) const {
    double p_xi = (field.p(i+1, j) - field.p(i-1, j)) / (2.0 * dxi_);
    double p_eta = (field.p(i, j+1) - field.p(i, j-1)) / (2.0 * deta_);
    return grid_.xi_y(i, j) * p_xi + grid_.eta_y(i, j) * p_eta;
}

// STEP 1: Advection-Diffusion
// SAME EQUATIONS as Cartesian, just with transformed derivatives!
//
// ∂u/∂t + u·∂u/∂x + v·∂u/∂y = ν·∇²u  (no pressure yet - that's step 2)
//
// STEP 1: Solve Advection-Diffusion with FULLY IMPLICIT SCHEME
//
// THEORY:
// For highly stretched grids, both advection AND diffusion must be implicit!
// 
// We solve:  u* = u + dt·[-∇·(u⊗u) + ν·∇²u]
//
// FULLY IMPLICIT (unconditionally stable):
// u* - dt·A(u*)·u* - dt·ν·∇²u* = u^n
// where A is the advection operator
//
// We linearize by evaluating velocities in advection operator at time n:
// [I - dt·A(u^n) - dt·ν·L]·u* = u^n
//
// This gives a LINEAR SYSTEM with both advection and diffusion implicit!
// Use first-order UPWIND for advection (stable, directional)
//
void CurvilinearSolver::computeIntermediateVelocities(CurvilinearField& field) {
    int nxi = grid_.nxi();
    int neta = grid_.neta();
    
    // STEP 1A: Compute RHS = u^n (just current velocity, no explicit terms!)
    // All advection and diffusion will be implicit
    std::vector<double> rhs_u(nxi * neta, 0.0);
    std::vector<double> rhs_v(nxi * neta, 0.0);
    
    // RHS is simply current velocity
    for (int j = 0; j < neta; j++) {
        for (int i = 0; i < nxi; i++) {
            int idx = j * nxi + i;
            rhs_u[idx] = field.u(i, j);
            rhs_v[idx] = field.v(i, j);
        }
    }
    
    // STEP 1B: Solve fully implicit system [I - dt·A - dt·ν·L]·u* = u^n
    // where A is upwind advection operator, L is diffusion operator
    const int max_iter = 100;  // More iterations needed for implicit advection
    const double tol = 1e-6;
    
    // Initialize u*, v* with current velocity as initial guess
    for (int j = 1; j < neta - 1; j++) {
        for (int i = 1; i < nxi - 1; i++) {
            field.u_star(i, j) = field.u(i, j);
            field.v_star(i, j) = field.v(i, j);
        }
    }
    
    // CRITICAL: Apply boundary conditions to u*, v* BEFORE iteration
    // Airfoil surface: no-slip
    for (int i = 0; i < nxi; i++) {
        field.u_star(i, 0) = 0.0;
        field.v_star(i, 0) = 0.0;
    }
    // Farfield: freestream
    for (int i = 0; i < nxi; i++) {
        field.u_star(i, neta-1) = 5.0;
        field.v_star(i, neta-1) = 0.0;
    }
    // Wake cut: extrapolate from interior
    for (int j = 0; j < neta; j++) {
        field.u_star(0, j) = field.u_star(1, j);
        field.v_star(0, j) = field.v_star(1, j);
        field.u_star(nxi-1, j) = field.u_star(nxi-2, j);
        field.v_star(nxi-1, j) = field.v_star(nxi-2, j);
    }
    
    // Jacobi iteration for implicit system
    // At each iteration: u*_new = RHS + dt·[A(u*_old) + ν·L(u*_old)]
    for (int iter = 0; iter < max_iter; iter++) {
        double max_change_u = 0.0;
        double max_change_v = 0.0;
        
        // Temporary storage for new values
        std::vector<double> u_new(nxi * neta);
        std::vector<double> v_new(nxi * neta);
        
        for (int j = 1; j < neta - 1; j++) {
            for (int i = 1; i < nxi - 1; i++) {
                int idx = j * nxi + i;
                
                // Get current velocity for advection operator
                double u_ij = field.u(i, j);
                double v_ij = field.v(i, j);
                
                // Transform to contravariant velocities in computational space
                // U = u·ξ_x + v·ξ_y (contravariant velocity in ξ direction)
                // V = u·η_x + v·η_y (contravariant velocity in η direction)
                double xi_x = grid_.xi_x(i, j);
                double xi_y = grid_.xi_y(i, j);
                double eta_x = grid_.eta_x(i, j);
                double eta_y = grid_.eta_y(i, j);
                
                double U_contra = u_ij * xi_x + v_ij * xi_y;
                double V_contra = u_ij * eta_x + v_ij * eta_y;
                
                // UPWIND ADVECTION for u*
                // If U > 0, flow moves in +ξ direction → use backward difference
                // If U < 0, flow moves in -ξ direction → use forward difference
                double u_star_adv_xi, v_star_adv_xi;
                if (U_contra > 0) {
                    u_star_adv_xi = (field.u_star(i,j) - field.u_star(i-1,j)) / dxi_;
                    v_star_adv_xi = (field.v_star(i,j) - field.v_star(i-1,j)) / dxi_;
                } else {
                    u_star_adv_xi = (field.u_star(i+1,j) - field.u_star(i,j)) / dxi_;
                    v_star_adv_xi = (field.v_star(i+1,j) - field.v_star(i,j)) / dxi_;
                }
                
                double u_star_adv_eta, v_star_adv_eta;
                if (V_contra > 0) {
                    u_star_adv_eta = (field.u_star(i,j) - field.u_star(i,j-1)) / deta_;
                    v_star_adv_eta = (field.v_star(i,j) - field.v_star(i,j-1)) / deta_;
                } else {
                    u_star_adv_eta = (field.u_star(i,j+1) - field.u_star(i,j)) / deta_;
                    v_star_adv_eta = (field.v_star(i,j+1) - field.v_star(i,j)) / deta_;
                }
                
                // Implicit advection contribution: -dt·(U·∂u*/∂ξ + V·∂u*/∂η)
                double adv_u_implicit = -dt_ * (U_contra * u_star_adv_xi + V_contra * u_star_adv_eta);
                double adv_v_implicit = -dt_ * (U_contra * v_star_adv_xi + V_contra * v_star_adv_eta);
                
                // Implicit diffusion contribution: dt·ν·∇²u*
                double diff_u_implicit = dt_ * nu_ * laplacian_u_star(field, i, j);
                double diff_v_implicit = dt_ * nu_ * laplacian_v_star(field, i, j);
                
                // New value: u*_new = RHS + dt·[advection + diffusion] evaluated at u*_old
                u_new[idx] = rhs_u[idx] + adv_u_implicit + diff_u_implicit;
                v_new[idx] = rhs_v[idx] + adv_v_implicit + diff_v_implicit;
                
                max_change_u = std::max(max_change_u, std::abs(u_new[idx] - field.u_star(i,j)));
                max_change_v = std::max(max_change_v, std::abs(v_new[idx] - field.v_star(i,j)));
            }
        }
        
        // Update u*, v* for interior points
        for (int j = 1; j < neta - 1; j++) {
            for (int i = 1; i < nxi - 1; i++) {
                int idx = j * nxi + i;
                field.u_star(i, j) = u_new[idx];
                field.v_star(i, j) = v_new[idx];
            }
        }
        
        // Update wake boundary with extrapolation
        for (int j = 1; j < neta - 1; j++) {
            field.u_star(0, j) = field.u_star(1, j);
            field.v_star(0, j) = field.v_star(1, j);
            field.u_star(nxi-1, j) = field.u_star(nxi-2, j);
            field.v_star(nxi-1, j) = field.v_star(nxi-2, j);
        }
        
        // Check convergence
        if (max_change_u < tol && max_change_v < tol) {
            break;
        }
    }
}

void CurvilinearSolver::applyBoundaryConditions(CurvilinearField& field) {
    int nxi = grid_.nxi();
    int neta = grid_.neta();
    
    // Airfoil surface (j=0): No-slip
    for (int i = 0; i < nxi; i++) {
        field.u(i, 0) = 0.0;
        field.v(i, 0) = 0.0;
        field.u_star(i, 0) = 0.0;
        field.v_star(i, 0) = 0.0;
    }
    
    // Farfield (j=neta-1): Freestream
    for (int i = 0; i < nxi; i++) {
        field.u(i, neta-1) = 5.0;  // U_inf
        field.v(i, neta-1) = 0.0;
        field.u_star(i, neta-1) = 5.0;
        field.v_star(i, neta-1) = 0.0;
    }
    
    // Wake cut (i=0 and i=nxi-1): Extrapolate from interior
    // Don't overwrite j=0 (airfoil) and j=neta-1 (farfield) which have their own BCs
    for (int j = 1; j < neta - 1; j++) {
        field.u(0, j) = field.u(1, j);
        field.v(0, j) = field.v(1, j);
        field.u(nxi-1, j) = field.u(nxi-2, j);
        field.v(nxi-1, j) = field.v(nxi-2, j);
        
        field.u_star(0, j) = field.u_star(1, j);
        field.v_star(0, j) = field.v_star(1, j);
        field.u_star(nxi-1, j) = field.u_star(nxi-2, j);
        field.v_star(nxi-1, j) = field.v_star(nxi-2, j);
    }
}

// STEP 2: Solve Pressure Poisson Equation
//
// THEORY:
// From incompressibility: ∇·u = 0
// After step 1 we have u*, which doesn't satisfy this.
// We need to find pressure p such that:
//   u^(n+1) = u* - dt·∇p  satisfies  ∇·u^(n+1) = 0
//
// This gives: ∇²p = (1/dt)·∇·u*
//
// IN CURVILINEAR COORDINATES (Thompson Eq. 5.18):
// The Laplacian ∇²p transforms to:
//   α·p_ξξ + 2β·p_ξη + γ·p_ηη + (metric derivatives) = RHS
//
// where RHS = (1/dt) · divergence of u* in transformed coords
//
void CurvilinearSolver::solvePressurePoisson(CurvilinearField& field) {
    int nxi = grid_.nxi();
    int neta = grid_.neta();
    
    const int max_iter = 200;
    const double omega = 1.5;  // SOR relaxation
    const double tolerance = 1e-6;
    
    // Iterative solver (SOR - same as Cartesian, but with transformed Laplacian)
    for (int iter = 0; iter < max_iter; iter++) {
        double max_change = 0.0;
        
        for (int j = 1; j < neta - 1; j++) {
            for (int i = 1; i < nxi - 1; i++) {
                double p_old = field.p(i, j);
                
                // Compute divergence of u* (RHS of Poisson equation)
                // ∇·u* = ∂u*/∂x + ∂v*/∂y
                // 
                // IN CURVILINEAR (Thompson Eq. 5.19):
                // Divergence transforms as:
                // ∇·u = (1/J)·[∂(J·U)/∂ξ + ∂(J·V)/∂η]
                // where U, V are contravariant velocity components
                //
                // For simplicity, use direct transformation:
                double u_star_xi = (field.u_star(i+1,j) - field.u_star(i-1,j)) / (2.0*dxi_);
                double u_star_eta = (field.u_star(i,j+1) - field.u_star(i,j-1)) / (2.0*deta_);
                double v_star_xi = (field.v_star(i+1,j) - field.v_star(i-1,j)) / (2.0*dxi_);
                double v_star_eta = (field.v_star(i,j+1) - field.v_star(i,j-1)) / (2.0*deta_);
                
                // Transform to physical space
                double du_star_dx = grid_.xi_x(i,j)*u_star_xi + grid_.eta_x(i,j)*u_star_eta;
                double dv_star_dy = grid_.xi_y(i,j)*v_star_xi + grid_.eta_y(i,j)*v_star_eta;
                
                double div_u_star = du_star_dx + dv_star_dy;
                
                // Safety check - if divergence too large, something's wrong
                if (std::abs(div_u_star) > 1e6) {
                    div_u_star = 0.0;  // Clamp to avoid blowup
                }
                
                double RHS = div_u_star / dt_;
                
                // Get neighbors
                double p_ip = field.p(i+1, j);
                double p_im = field.p(i-1, j);
                double p_jp = field.p(i, j+1);
                double p_jm = field.p(i, j-1);
                
                // Corner points for mixed derivative
                double p_ipjp = field.p(i+1, j+1);
                double p_imjp = field.p(i-1, j+1);
                double p_ipjm = field.p(i+1, j-1);
                double p_imjm = field.p(i-1, j-1);
                
                // Compute metric coefficients at this point
                double xi_x = grid_.xi_x(i, j);
                double xi_y = grid_.xi_y(i, j);
                double eta_x = grid_.eta_x(i, j);
                double eta_y = grid_.eta_y(i, j);
                
                double alpha = xi_x*xi_x + xi_y*xi_y;
                double beta = xi_x*eta_x + xi_y*eta_y;
                double gamma = eta_x*eta_x + eta_y*eta_y;
                
                // Discretized Laplacian operator (Thompson Eq. 5.20)
                // α·p_ξξ + 2β·p_ξη + γ·p_ηη = RHS
                //
                // p_ξξ ≈ (p[i+1] - 2p[i] + p[i-1])/dξ²
                // p_ηη ≈ (p[j+1] - 2p[j] + p[j-1])/dη²
                // p_ξη ≈ (p[i+1,j+1] - p[i+1,j-1] - p[i-1,j+1] + p[i-1,j-1])/(4dξdη)
                
                double coef_center = 2.0 * (alpha/(dxi_*dxi_) + gamma/(deta_*deta_));
                double coef_xi = alpha / (dxi_*dxi_);
                double coef_eta = gamma / (deta_*deta_);
                double coef_xieta = beta / (2.0*dxi_*deta_);
                
                // Solve for new pressure (rearranging the discrete Laplacian)
                double p_new = (coef_xi * (p_ip + p_im)
                              + coef_eta * (p_jp + p_jm)
                              + coef_xieta * (p_ipjp - p_ipjm - p_imjp + p_imjm)
                              - RHS) / coef_center;
                
                // SOR update
                field.p(i, j) = (1.0 - omega) * p_old + omega * p_new;
                
                max_change = std::max(max_change, std::abs(field.p(i,j) - p_old));
            }
        }
        
        // Pressure boundary conditions
        // Airfoil surface: ∂p/∂η = 0 (zero normal gradient)
        for (int i = 0; i < nxi; i++) {
            field.p(i, 0) = field.p(i, 1);
        }
        
        // Farfield: ∂p/∂η = 0
        for (int i = 0; i < nxi; i++) {
            field.p(i, neta-1) = field.p(i, neta-2);
        }
        
        // Wake cut: extrapolate
        for (int j = 0; j < neta; j++) {
            field.p(0, j) = field.p(1, j);
            field.p(nxi-1, j) = field.p(nxi-2, j);
        }
        
        // Check convergence
        if (max_change < tolerance) {
            break;
        }
    }
    
    // Remove pressure mean (arbitrary constant)
    double p_sum = 0.0;
    int count = 0;
    for (int j = 1; j < neta-1; j++) {
        for (int i = 1; i < nxi-1; i++) {
            p_sum += field.p(i, j);
            count++;
        }
    }
    double p_mean = p_sum / count;
    
    for (int j = 0; j < neta; j++) {
        for (int i = 0; i < nxi; i++) {
            field.p(i, j) -= p_mean;
        }
    }
}

// STEP 3: Velocity Correction
//
// THEORY:
// u^(n+1) = u* - dt·∇p
//
// IN CURVILINEAR:
// ∇p is transformed using metrics (same as du_dx, etc.)
//
void CurvilinearSolver::correctVelocities(CurvilinearField& field) {
    int nxi = grid_.nxi();
    int neta = grid_.neta();
    
    for (int j = 1; j < neta - 1; j++) {
        for (int i = 1; i < nxi - 1; i++) {
            // Skip airfoil (will be overwritten by BC)
            if (j == 0) continue;
            
            // Compute pressure gradient in physical space
            // Uses transformed derivatives with metrics
            double dpdx = dp_dx(field, i, j);
            double dpdy = dp_dy(field, i, j);
            
            // Correct velocities
            field.u(i, j) = field.u_star(i, j) - dt_ * dpdx;
            field.v(i, j) = field.v_star(i, j) - dt_ * dpdy;
        }
    }
}

void CurvilinearSolver::step(CurvilinearField& field) {
    applyBoundaryConditions(field);
    computeIntermediateVelocities(field);
    solvePressurePoisson(field);
    correctVelocities(field);
    applyBoundaryConditions(field);
}
