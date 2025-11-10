#pragma once
#include <vector>
#include <string>
#include <cmath>

// Body-fitted C-grid generator for airfoils using elliptic PDE method
// Based on:
//   - Sorenson (1980) NASA TM 81198 - Inverse Poisson formulation
//   - Thompson et al. (1985) Ch. 5 - Elliptic generation theory
//   - Steger & Sorenson (1979) - Boundary clustering

class MeshGenerator {
public:
    // Grid dimensions in computational space (ξ, η)
    MeshGenerator(int nxi, int neta);
    
    // Generate C-grid around NACA 4-digit airfoil
    // chord: airfoil chord length (m)
    // alpha_deg: angle of attack (degrees)
    // farfield_radius: distance to outer boundary (typically 15-20 chords)
    void generateCGrid(double chord, double alpha_deg, double farfield_radius);
    
    // Generate O-grid around NACA 4-digit airfoil (better topology - no wake cut!)
    // Same parameters as C-grid but wraps smoothly around airfoil
    void generateOGrid(double chord, double alpha_deg, double farfield_radius);
    
    // Get physical coordinates at grid point (i,j)
    double x(int i, int j) const { return x_[j * nxi_ + i]; }
    double y(int i, int j) const { return y_[j * nxi_ + i]; }
    
    // Get grid metrics (needed for transformed Navier-Stokes)
    double J(int i, int j) const { return J_[j * nxi_ + i]; }       // Jacobian
    double xi_x(int i, int j) const { return xi_x_[j * nxi_ + i]; } // ∂ξ/∂x
    double xi_y(int i, int j) const { return xi_y_[j * nxi_ + i]; } // ∂ξ/∂y
    double eta_x(int i, int j) const { return eta_x_[j * nxi_ + i]; } // ∂η/∂x
    double eta_y(int i, int j) const { return eta_y_[j * nxi_ + i]; } // ∂η/∂y
    
    // Grid dimensions
    int nxi() const { return nxi_; }
    int neta() const { return neta_; }
    
    // Export grid for visualization
    void exportGrid(const std::string& filename) const;
    void exportMetrics(const std::string& filename) const;

private:
    int nxi_;   // Grid points in ξ direction (wraps around airfoil)
    int neta_;  // Grid points in η direction (normal to airfoil)
    
    // Physical coordinates (x,y) at each (ξ,η) point
    std::vector<double> x_, y_;
    
    // Grid metrics at each point
    std::vector<double> J_;      // Jacobian = x_ξ*y_η - x_η*y_ξ
    std::vector<double> xi_x_;   // ∂ξ/∂x = y_η/J
    std::vector<double> xi_y_;   // ∂ξ/∂y = -x_η/J
    std::vector<double> eta_x_;  // ∂η/∂x = -y_ξ/J
    std::vector<double> eta_y_;  // ∂η/∂y = x_ξ/J
    
    // NACA 4-digit airfoil geometry
    // From: Abbott & Von Doenhoff (1959) "Theory of Wing Sections"
    // Returns (x,y) on airfoil surface for parameter t ∈ [0,1]
    // t=0: trailing edge (upper), t=0.5: leading edge, t=1: trailing edge (lower)
    void nacaAirfoil(double t, double chord, double& x, double& y) const;
    
    // Step 1: Algebraic initialization (Sorenson Eq. 3)
    // Simple transfinite interpolation between boundaries
    // Gets topology right but not smooth/orthogonal
    void initializeAlgebraic(double chord, double alpha_deg, double farfield_radius);
    
    // Step 2: Elliptic smoothing (Thompson Ch. 5, Sorenson Eq. 15-16)
    // Solves inverse Poisson equations:
    //   α·x_ξξ - 2β·x_ξη + γ·x_ηη = -J²(P·x_ξ + Q·x_η)
    // where P,Q are control functions (start with P=Q=0 for pure Laplace)
    void smoothElliptic(int max_iterations, double tolerance);
    
    // Step 3: Compute metrics from final grid (used in CFD solver)
    // Central differences for derivatives
    void computeMetrics();
    
    // Helper: Control functions for orthogonality (Thomas & Middlecoff 1980)
    // Start with P=Q=0, add later for better quality
    double controlP(int i, int j) const { return 0.0; }
    double controlQ(int i, int j) const { return 0.0; }
};
