#include "mesh_generator.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

MeshGenerator::MeshGenerator(int nxi, int neta) 
    : nxi_(nxi), neta_(neta) {
    // Allocate storage for grid points and metrics
    int total_points = nxi_ * neta_;
    x_.resize(total_points);
    y_.resize(total_points);
    J_.resize(total_points);
    xi_x_.resize(total_points);
    xi_y_.resize(total_points);
    eta_x_.resize(total_points);
    eta_y_.resize(total_points);
}

// NACA 4-digit airfoil equations
// Reference: Abbott & Von Doenhoff (1959) "Theory of Wing Sections" p. 112-113
// NACA 2412: 2% max camber at 40% chord, 12% thickness
void MeshGenerator::nacaAirfoil(double t, double chord, double& xout, double& yout) const {
    // For NACA 2412 (hardcoded for now, can parameterize later)
    const double m = 0.02;  // Max camber (2%)
    const double p = 0.4;   // Position of max camber (40% chord)
    const double thick = 0.12; // Max thickness (12%)
    
    // Parameter t ∈ [0,1] wraps around airfoil:
    //   t=0.0   → Trailing edge (upper surface)
    //   t=0.25  → Mid-chord upper
    //   t=0.5   → Leading edge
    //   t=0.75  → Mid-chord lower  
    //   t=1.0   → Trailing edge (lower surface)
    
    double x_c; // Chordwise position [0,1]
    bool upper_surface;
    
    if (t <= 0.5) {
        // Upper surface: parametrize from TE to LE
        x_c = 1.0 - 2.0 * t;  // t=0→x=1, t=0.5→x=0
        upper_surface = true;
    } else {
        // Lower surface: parametrize from LE to TE
        x_c = 2.0 * (t - 0.5); // t=0.5→x=0, t=1→x=1
        upper_surface = false;
    }
    
    // Thickness distribution (NACA 4-digit, Eq. 6.4 from Abbott & Von Doenhoff)
    // Coefficients give max thickness at ~30% chord
    double yt = 5.0 * thick * (
        0.2969 * std::sqrt(x_c) 
        - 0.1260 * x_c 
        - 0.3516 * x_c * x_c 
        + 0.2843 * x_c * x_c * x_c 
        - 0.1015 * x_c * x_c * x_c * x_c
    );
    
    // Camber line (NACA 4-digit, Eq. 6.1-6.2)
    double yc, dyc_dx;
    if (x_c < p) {
        yc = m / (p * p) * (2.0 * p * x_c - x_c * x_c);
        dyc_dx = 2.0 * m / (p * p) * (p - x_c);
    } else {
        yc = m / ((1.0 - p) * (1.0 - p)) * (1.0 - 2.0 * p + 2.0 * p * x_c - x_c * x_c);
        dyc_dx = 2.0 * m / ((1.0 - p) * (1.0 - p)) * (p - x_c);
    }
    
    // Angle of camber line
    double theta = std::atan(dyc_dx);
    
    // Add/subtract thickness perpendicular to camber line (Eq. 6.3)
    if (upper_surface) {
        xout = chord * (x_c - yt * std::sin(theta));
        yout = chord * (yc + yt * std::cos(theta));
    } else {
        xout = chord * (x_c + yt * std::sin(theta));
        yout = chord * (yc - yt * std::cos(theta));
    }
}

void MeshGenerator::generateCGrid(double chord, double alpha_deg, double farfield_radius) {
    std::cout << "Generating C-grid around NACA 2412 airfoil...\n";
    std::cout << "  Grid size: " << nxi_ << " x " << neta_ << " = " << nxi_*neta_ << " points\n";
    std::cout << "  Chord: " << chord << " m\n";
    std::cout << "  Angle of attack: " << alpha_deg << " deg\n";
    std::cout << "  Farfield radius: " << farfield_radius << " m\n\n";
    
    // Step 1: Initialize with algebraic transfinite interpolation
    // This is our "first guess" - gets the C-topology right
    std::cout << "Step 1/3: Algebraic initialization...\n";
    initializeAlgebraic(chord, alpha_deg, farfield_radius);
    
    // Step 2: Smooth with elliptic solver
    // This makes grid orthogonal at boundaries and removes kinks
    std::cout << "Step 2/3: Elliptic smoothing...\n";
    smoothElliptic(10000, 5e-8);  // Back to 10k iterations
    
    // Step 3: Compute metrics for CFD solver
    std::cout << "Step 3/3: Computing grid metrics...\n";
    computeMetrics();
    
    std::cout << "Grid generation complete!\n\n";
}

// STEP 1: Algebraic Initialization
// Reference: Sorenson (1980) Eq. 3, Thompson (1985) p. 147
// Why: Gets C-topology correct, serves as initial guess for elliptic solver
void MeshGenerator::initializeAlgebraic(double chord, double alpha_deg, double farfield_radius) {
    const double alpha_rad = alpha_deg * M_PI / 180.0;
    
    // Grid stretching parameter (controls wall clustering)
    // β=0: uniform spacing, β>0: clustering near airfoil (η=0)
    // Using tanh stretching prevents extreme cell aspect ratios
    const double beta = 0.0;  // DISABLED - linear spacing works better for this grid
    
    std::cout << "  Using " << (beta > 0.01 ? "tanh stretching" : "uniform spacing") 
              << " in η-direction\n";
    
    // Define boundaries:
    // η=0: Airfoil surface (parametrized by ξ)
    // η=1: Far-field boundary (large C-shape)
    
    for (int j = 0; j < neta_; j++) {
        double eta_uniform = j / double(neta_ - 1);  // Uniform η ∈ [0,1]
        
        // Apply hyperbolic tangent stretching (Thompson & Mastin 1985)
        // This clusters points near η=0 (airfoil) while maintaining smoothness
        double eta;
        if (beta > 0.01) {
            eta = std::tanh(beta * eta_uniform) / std::tanh(beta);
        } else {
            eta = eta_uniform;  // Fall back to uniform if β≈0
        }
        
        for (int i = 0; i < nxi_; i++) {
            double xi = i / double(nxi_ - 1);  // ξ ∈ [0,1]
            
            // Boundary at η=0 (airfoil surface)
            double x_inner, y_inner;
            nacaAirfoil(xi, chord, x_inner, y_inner);
            
            // Boundary at η=1 (far-field C-shape)
            // The "C" wraps around airfoil - use a simple circular arc
            // Reference: Thompson (1985) Fig. 7.3, p. 236 - "C-Grid for Airfoil"
            double x_outer, y_outer;
            
            // Map ξ ∈ [0,1] to angle θ around circle
            // ξ=0 & ξ=1: Downstream wake (θ=±π)
            // ξ=0.5: Leading edge side (θ=0)
            double theta;
            if (xi <= 0.5) {
                // Upper half: ξ=0 → θ=π, ξ=0.5 → θ=0
                theta = M_PI * (1.0 - 2.0 * xi);
            } else {
                // Lower half: ξ=0.5 → θ=0, ξ=1 → θ=-π
                theta = -M_PI * (2.0 * xi - 1.0);
            }
            
            // Circular boundary centered at mid-chord
            double x_center = 0.5 * chord;
            double y_center = 0.0;
            x_outer = x_center + farfield_radius * std::cos(theta);
            y_outer = y_center + farfield_radius * std::sin(theta);
            
            // Transfinite interpolation (linear blending, Sorenson Eq. 3)
            // x(ξ,η) = (1-η)·x_inner(ξ) + η·x_outer(ξ)
            double x_blend = (1.0 - eta) * x_inner + eta * x_outer;
            double y_blend = (1.0 - eta) * y_inner + eta * y_outer;
            
            // Apply angle of attack rotation
            // Rotate grid so freestream is horizontal
            double x_rot = x_blend * std::cos(alpha_rad) + y_blend * std::sin(alpha_rad);
            double y_rot = -x_blend * std::sin(alpha_rad) + y_blend * std::cos(alpha_rad);
            
            x_[j * nxi_ + i] = x_rot;
            y_[j * nxi_ + i] = y_rot;
        }
    }
    
    std::cout << "  Algebraic grid initialized with C-topology\n";
}

void MeshGenerator::exportGrid(const std::string& filename) const {
    std::ofstream file(filename);
    file << "x,y,i,j\n";
    for (int j = 0; j < neta_; j++) {
        for (int i = 0; i < nxi_; i++) {
            file << x(i,j) << "," << y(i,j) << "," << i << "," << j << "\n";
        }
    }
    std::cout << "Grid exported to " << filename << "\n";
}

void MeshGenerator::exportMetrics(const std::string& filename) const {
    std::ofstream file(filename);
    file << "i,j,J,xi_x,xi_y,eta_x,eta_y\n";
    for (int j = 0; j < neta_; j++) {
        for (int i = 0; i < nxi_; i++) {
            file << i << "," << j << "," 
                 << J(i,j) << "," 
                 << xi_x(i,j) << "," << xi_y(i,j) << "," 
                 << eta_x(i,j) << "," << eta_y(i,j) << "\n";
        }
    }
    std::cout << "Metrics exported to " << filename << "\n";
}

// STEP 2: Elliptic Smoothing with Inverse Poisson Equations
// Reference: Sorenson (1980) NASA TM 81198, Eq. 15-16
//            Thompson et al. (1985) Ch. 5, p. 162-165
//
// WHY THIS WORKS:
// We want grid lines to satisfy Laplace equation (like soap film):
//   ∇²ξ = 0  and  ∇²η = 0  in physical space (x,y)
//
// But we don't know (x,y) yet! We only know computational space (ξ,η).
// So we INVERT the problem using chain rule (Thompson p. 163):
//
// Original:  ∂²ξ/∂x² + ∂²ξ/∂y² = 0   (easy if you know x,y)
// Inverted:  α·∂²x/∂ξ² - 2β·∂²x/∂ξ∂η + γ·∂²x/∂η² = 0  (easy if you know ξ,η)
//
// where α, β, γ are metric coefficients computed from current grid.
// Same equations for y. This is genius because our computational grid
// is a perfect rectangle!
//
// ALGORITHM (Sorenson Section 3, Thompson p. 165):
// 1. Start with algebraic grid (already done)
// 2. Compute metrics α, β, γ from current x(ξ,η), y(ξ,η)
// 3. Solve for new x, y using SOR (Successive Over-Relaxation)
// 4. Repeat until convergence (grid stops changing)
//
void MeshGenerator::smoothElliptic(int max_iterations, double tolerance) {
    const double dxi = 1.0 / (nxi_ - 1);
    const double deta = 1.0 / (neta_ - 1);
    const double omega = 1.0;  // Start with Gauss-Seidel (ω=1.0)
                                // Will increase to 1.5 once stable
    
    std::cout << "  Elliptic solver settings:\n";
    std::cout << "    Max iterations: " << max_iterations << "\n";
    std::cout << "    Tolerance: " << tolerance << "\n";
    std::cout << "    SOR omega: " << omega << "\n\n";
    
    // Precompute constants for efficiency
    const double dxi2 = dxi * dxi;
    const double deta2 = deta * deta;
    const double dxi_deta = dxi * deta;
    
    for (int iter = 0; iter < max_iterations; iter++) {
        double max_update = 0.0;
        int nan_count = 0;
        
        // Sweep interior points (boundaries are fixed)
        // Use "point SOR" - update each point using most recent neighbor values
        // Reference: Sorenson (1980) Section 3, Thompson (1985) p. 165-166
        for (int j = 1; j < neta_ - 1; j++) {
            for (int i = 1; i < nxi_ - 1; i++) {
                int idx = j * nxi_ + i;
                
                // Store old values
                double x_old = x_[idx];
                double y_old = y_[idx];
                
                // Get neighbor values (some may already be updated in this iteration)
                double x_ip = x_[idx+1];      // i+1
                double x_im = x_[idx-1];      // i-1
                double x_jp = x_[idx+nxi_];   // j+1
                double x_jm = x_[idx-nxi_];   // j-1
                
                double y_ip = y_[idx+1];
                double y_im = y_[idx-1];
                double y_jp = y_[idx+nxi_];
                double y_jm = y_[idx-nxi_];
                
                // Compute metrics from current grid
                // (Using old value at center, newest values at neighbors)
                double x_xi = (x_ip - x_im) / (2.0 * dxi);
                double x_eta = (x_jp - x_jm) / (2.0 * deta);
                double y_xi = (y_ip - y_im) / (2.0 * dxi);
                double y_eta = (y_jp - y_jm) / (2.0 * deta);
                
                // Metric coefficients
                double alpha = x_eta*x_eta + y_eta*y_eta;
                double beta = x_xi*x_eta + y_xi*y_eta;
                double gamma = x_xi*x_xi + y_xi*y_xi;
                
                // Solve for new x, y from discretized Poisson equation
                // α·x_ξξ - 2β·x_ξη + γ·x_ηη = 0
                //
                // Discretization (Sorenson Eq. 17-18):
                // x[i][j] = [α·(x[i+1][j] + x[i-1][j])/dξ² 
                //          + γ·(x[i][j+1] + x[i][j-1])/dη²
                //          - β·(x[i+1][j+1] - x[i+1][j-1] - x[i-1][j+1] + x[i-1][j-1])/(2·dξ·dη)]
                //         / [2(α/dξ² + γ/dη²)]
                
                double denominator = 2.0 * (alpha/dxi2 + gamma/deta2);
                
                if (std::abs(denominator) < 1e-10) {
                    // Skip if denominator too small (shouldn't happen with valid grid)
                    continue;
                }
                
                // Mixed derivative term
                double x_xieta_contrib = beta * (x_[idx+nxi_+1] - x_[idx+nxi_-1] 
                                                - x_[idx-nxi_+1] + x_[idx-nxi_-1]) 
                                        / (2.0 * dxi_deta);
                double y_xieta_contrib = beta * (y_[idx+nxi_+1] - y_[idx+nxi_-1] 
                                                - y_[idx-nxi_+1] + y_[idx-nxi_-1]) 
                                        / (2.0 * dxi_deta);
                
                // New values from Laplace equation
                double x_new = (alpha * (x_ip + x_im) / dxi2
                              + gamma * (x_jp + x_jm) / deta2
                              - x_xieta_contrib) / denominator;
                              
                double y_new = (alpha * (y_ip + y_im) / dxi2
                              + gamma * (y_jp + y_jm) / deta2
                              - y_xieta_contrib) / denominator;
                
                // Apply SOR relaxation
                x_[idx] = (1.0 - omega) * x_old + omega * x_new;
                y_[idx] = (1.0 - omega) * y_old + omega * y_new;
                
                // Check for NaN/Inf
                if (std::isnan(x_[idx]) || std::isnan(y_[idx]) || 
                    std::isinf(x_[idx]) || std::isinf(y_[idx])) {
                    nan_count++;
                    x_[idx] = x_old;
                    y_[idx] = y_old;
                }
                
                // Track maximum change
                double update = std::sqrt((x_[idx] - x_old)*(x_[idx] - x_old) + 
                                         (y_[idx] - y_old)*(y_[idx] - y_old));
                max_update = std::max(max_update, update);
            }
        }
        
        // Print progress
        if (iter % 100 == 0 || iter < 10) {
            std::cout << "    Iteration " << iter 
                      << ": max update = " << max_update;
            if (nan_count > 0) {
                std::cout << " [" << nan_count << " NaNs]";
            }
            std::cout << "\n";
        }
        
        // Check for fatal errors
        if (nan_count > 10) {
            std::cout << "\n  ⚠ ERROR: NaN values detected - grid generation failed\n\n";
            return;
        }
        
        // Check convergence
        if (max_update < tolerance) {
            std::cout << "\n  ✓ Converged in " << iter + 1 << " iterations!\n";
            std::cout << "    Final max update: " << max_update << "\n\n";
            return;
        }
    }
    
    std::cout << "\n  ⚠ Did not fully converge in " << max_iterations << " iterations\n";
    std::cout << "    Grid may still be usable - check Jacobian\n\n";
}

// STEP 3: Compute Grid Metrics
// Reference: Thompson et al. (1985) Eq. 5.9-5.11, p. 164
//
// WHY WE NEED THESE:
// To transform Navier-Stokes from (x,y) to (ξ,η), we need:
//   ∂/∂x = ξₓ·∂/∂ξ + ηₓ·∂/∂η
//   ∂/∂y = ξᵧ·∂/∂ξ + ηᵧ·∂/∂η
//
// These ξₓ, ξᵧ, ηₓ, ηᵧ are called "inverse metrics" or "transformation metrics"
//
void MeshGenerator::computeMetrics() {
    const double dxi = 1.0 / (nxi_ - 1);
    const double deta = 1.0 / (neta_ - 1);
    
    double min_jacobian = 1e10;
    double max_jacobian = -1e10;
    int negative_count = 0;
    
    // Interior points (use central differences)
    for (int j = 1; j < neta_ - 1; j++) {
        for (int i = 1; i < nxi_ - 1; i++) {
            int idx = j * nxi_ + i;
            
            // Compute forward transformation metrics (Thompson Eq. 5.9)
            double x_xi = (x_[idx+1] - x_[idx-1]) / (2.0 * dxi);
            double x_eta = (x_[idx+nxi_] - x_[idx-nxi_]) / (2.0 * deta);
            double y_xi = (y_[idx+1] - y_[idx-1]) / (2.0 * dxi);
            double y_eta = (y_[idx+nxi_] - y_[idx-nxi_]) / (2.0 * deta);
            
            // Jacobian (Thompson Eq. 5.10)
            // J = det of transformation matrix = |∂(x,y)/∂(ξ,η)|
            // MUST be positive everywhere for valid grid!
            J_[idx] = x_xi * y_eta - x_eta * y_xi;
            
            if (J_[idx] <= 0.0) negative_count++;
            min_jacobian = std::min(min_jacobian, J_[idx]);
            max_jacobian = std::max(max_jacobian, J_[idx]);
            
            // Inverse metrics (Thompson Eq. 5.11)
            // These come from inverting the transformation Jacobian matrix
            xi_x_[idx] = y_eta / J_[idx];
            xi_y_[idx] = -x_eta / J_[idx];
            eta_x_[idx] = -y_xi / J_[idx];
            eta_y_[idx] = x_xi / J_[idx];
        }
    }
    
    // Boundary points (use one-sided differences)
    for (int i = 0; i < nxi_; i++) {
        // j=0 (airfoil surface)
        {
            int idx = 0 * nxi_ + i;
            double x_xi = i == 0 ? (x_[idx+1] - x_[idx]) / dxi :
                         i == nxi_-1 ? (x_[idx] - x_[idx-1]) / dxi :
                         (x_[idx+1] - x_[idx-1]) / (2.0 * dxi);
            double x_eta = (x_[idx+nxi_] - x_[idx]) / deta;
            double y_xi = i == 0 ? (y_[idx+1] - y_[idx]) / dxi :
                         i == nxi_-1 ? (y_[idx] - y_[idx-1]) / dxi :
                         (y_[idx+1] - y_[idx-1]) / (2.0 * dxi);
            double y_eta = (y_[idx+nxi_] - y_[idx]) / deta;
            
            J_[idx] = x_xi * y_eta - x_eta * y_xi;
            if (J_[idx] <= 0.0) negative_count++;
            min_jacobian = std::min(min_jacobian, J_[idx]);
            max_jacobian = std::max(max_jacobian, J_[idx]);
            
            xi_x_[idx] = y_eta / J_[idx];
            xi_y_[idx] = -x_eta / J_[idx];
            eta_x_[idx] = -y_xi / J_[idx];
            eta_y_[idx] = x_xi / J_[idx];
        }
        
        // j=neta-1 (farfield)
        {
            int idx = (neta_-1) * nxi_ + i;
            double x_xi = i == 0 ? (x_[idx+1] - x_[idx]) / dxi :
                         i == nxi_-1 ? (x_[idx] - x_[idx-1]) / dxi :
                         (x_[idx+1] - x_[idx-1]) / (2.0 * dxi);
            double x_eta = (x_[idx] - x_[idx-nxi_]) / deta;
            double y_xi = i == 0 ? (y_[idx+1] - y_[idx]) / dxi :
                         i == nxi_-1 ? (y_[idx] - y_[idx-1]) / dxi :
                         (y_[idx+1] - y_[idx-1]) / (2.0 * dxi);
            double y_eta = (y_[idx] - y_[idx-nxi_]) / deta;
            
            J_[idx] = x_xi * y_eta - x_eta * y_xi;
            if (J_[idx] <= 0.0) negative_count++;
            min_jacobian = std::min(min_jacobian, J_[idx]);
            max_jacobian = std::max(max_jacobian, J_[idx]);
            
            xi_x_[idx] = y_eta / J_[idx];
            xi_y_[idx] = -x_eta / J_[idx];
            eta_x_[idx] = -y_xi / J_[idx];
            eta_y_[idx] = x_xi / J_[idx];
        }
    }
    
    // i-boundaries (wake cut edges)
    for (int j = 1; j < neta_ - 1; j++) {
        // i=0 (left edge of wake cut)
        {
            int idx = j * nxi_ + 0;
            double x_xi = (x_[idx+1] - x_[idx]) / dxi;  // Forward diff
            double x_eta = (x_[idx+nxi_] - x_[idx-nxi_]) / (2.0 * deta);  // Central
            double y_xi = (y_[idx+1] - y_[idx]) / dxi;
            double y_eta = (y_[idx+nxi_] - y_[idx-nxi_]) / (2.0 * deta);
            
            J_[idx] = x_xi * y_eta - x_eta * y_xi;
            if (J_[idx] <= 0.0) negative_count++;
            min_jacobian = std::min(min_jacobian, J_[idx]);
            max_jacobian = std::max(max_jacobian, J_[idx]);
            
            xi_x_[idx] = y_eta / J_[idx];
            xi_y_[idx] = -x_eta / J_[idx];
            eta_x_[idx] = -y_xi / J_[idx];
            eta_y_[idx] = x_xi / J_[idx];
        }
        
        // i=nxi-1 (right edge of wake cut)
        {
            int idx = j * nxi_ + (nxi_ - 1);
            double x_xi = (x_[idx] - x_[idx-1]) / dxi;  // Backward diff
            double x_eta = (x_[idx+nxi_] - x_[idx-nxi_]) / (2.0 * deta);  // Central
            double y_xi = (y_[idx] - y_[idx-1]) / dxi;
            double y_eta = (y_[idx+nxi_] - y_[idx-nxi_]) / (2.0 * deta);
            
            J_[idx] = x_xi * y_eta - x_eta * y_xi;
            if (J_[idx] <= 0.0) negative_count++;
            min_jacobian = std::min(min_jacobian, J_[idx]);
            max_jacobian = std::max(max_jacobian, J_[idx]);
            
            xi_x_[idx] = y_eta / J_[idx];
            xi_y_[idx] = -x_eta / J_[idx];
            eta_x_[idx] = -y_xi / J_[idx];
            eta_y_[idx] = x_xi / J_[idx];
        }
    }
    
    // Validate grid quality (Sorenson Section 4.1)
    std::cout << "  Jacobian range: [" << min_jacobian << ", " << max_jacobian << "]\n";
    
    if (negative_count > 0) {
        std::cout << "  ⚠ WARNING: " << negative_count << " points have negative Jacobian!\n";
        std::cout << "    Grid has folded or crossed lines.\n";
        std::cout << "    Increase farfield radius or adjust clustering.\n";
    } else {
        std::cout << "  ✓ All Jacobians positive - grid is valid!\n";
    }
    
    // Check grid aspect ratio (should be O(1) for good quality)
    double aspect_ratio = max_jacobian / min_jacobian;
    std::cout << "  Grid aspect ratio: " << aspect_ratio << "\n";
    
    if (aspect_ratio > 100.0) {
        std::cout << "  ⚠ High aspect ratio - consider adjusting clustering\n";
    }
    
    std::cout << "\n";
}
