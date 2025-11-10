#pragma once
#include "curvilinear_grid.hpp"
#include "curvilinear_field.hpp"

// Navier-Stokes solver on body-fitted curvilinear grid
// 
// KEY DIFFERENCE FROM CARTESIAN:
// Instead of ∂u/∂x, we compute: ∂u/∂x = ξ_x·∂u/∂ξ + η_x·∂u/∂η
// 
// WHY THIS WORKS:
// The Navier-Stokes equations are the SAME physical laws.
// We're just using a different coordinate system (ξ,η) to describe them.
// 
// TRANSFORMATION (Thompson Ch. 5):
// ∂/∂x = ξ_x·∂/∂ξ + η_x·∂/∂η
// ∂/∂y = ξ_y·∂/∂ξ + η_y·∂/∂η
//
// LAPLACIAN becomes (Thompson Eq. 5.15):
// ∇²u = α·u_ξξ + 2β·u_ξη + γ·u_ηη + first derivative terms
// where α, β, γ are metric coefficients (same as grid generation!)
//
// Reference: Thompson et al. (1985) Ch. 5, Sections 5.4-5.5

class CurvilinearSolver {
public:
    CurvilinearSolver(const CurvilinearGrid& grid, double Re, double dt);
    
    // Take one timestep using fractional-step method
    // Same algorithm as Cartesian, but with transformed derivatives
    void step(CurvilinearField& field);
    
    // Step 1: Advection-diffusion (get u*, v*)
    void computeIntermediateVelocities(CurvilinearField& field);
    
    // Step 2: Pressure Poisson equation
    void solvePressurePoisson(CurvilinearField& field);
    
    // Step 3: Velocity correction
    void correctVelocities(CurvilinearField& field);
    
    // Apply boundary conditions on curvilinear grid
    void applyBoundaryConditions(CurvilinearField& field);
    
private:
    const CurvilinearGrid& grid_;
    double Re_;     // Reynolds number
    double nu_;     // Kinematic viscosity = 1/Re
    double dt_;     // Timestep
    
    double dxi_;    // Computational grid spacing (uniform!)
    double deta_;
    
    // Compute transformed derivative ∂u/∂x at point (i,j)
    // Uses metrics and finite differences in (ξ,η) space
    double du_dx(const CurvilinearField& field, int i, int j) const;
    double du_dy(const CurvilinearField& field, int i, int j) const;
    double dv_dx(const CurvilinearField& field, int i, int j) const;
    double dv_dy(const CurvilinearField& field, int i, int j) const;
    
    // Compute transformed Laplacian ∇²u
    double laplacian_u(const CurvilinearField& field, int i, int j) const;
    double laplacian_v(const CurvilinearField& field, int i, int j) const;
    
    // Laplacian of u_star, v_star (for implicit diffusion)
    double laplacian_u_star(const CurvilinearField& field, int i, int j) const;
    double laplacian_v_star(const CurvilinearField& field, int i, int j) const;
    
    // Pressure gradient in physical space
    double dp_dx(const CurvilinearField& field, int i, int j) const;
    double dp_dy(const CurvilinearField& field, int i, int j) const;
};
