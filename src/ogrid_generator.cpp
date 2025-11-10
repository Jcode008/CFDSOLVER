#include "mesh_generator.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// O-GRID GENERATOR
// 
// TOPOLOGY:
// - ξ: wraps around airfoil (0 → 2π, periodic)
// - η: radial direction from surface to farfield (0 → 1)
// 
// ADVANTAGES over C-grid:
// - No wake cut discontinuity
// - More uniform aspect ratio
// - Better for elliptic smoothing
//
// CHALLENGE:
// - Sharp trailing edge needs special treatment
//   Solution: Use small "sting" extension or slight bluntness
//

void MeshGenerator::generateOGrid(double chord, double alpha_deg, double farfield_radius) {
    std::cout << "Generating O-grid around NACA 2412 airfoil...\n";
    std::cout << "  Grid size: " << nxi_ << " x " << neta_ << " = " << nxi_*neta_ << " points\n";
    std::cout << "  Chord: " << chord << " m\n";
    std::cout << "  Angle of attack: " << alpha_deg << " deg\n";
    std::cout << "  Farfield radius: " << farfield_radius << " m\n\n";
    
    // NACA 2412 parameters
    const double m = 0.02;     // max camber (2%)
    const double p = 0.4;      // location of max camber (40% chord)
    const double thick = 0.12; // thickness (12%)
    
    // Store angle of attack for rotation
    double alpha = alpha_deg * M_PI / 180.0;
    
    // STEP 1: Generate airfoil surface points
    std::cout << "Step 1/3: Generating airfoil surface...\n";
    std::vector<double> x_airfoil(nxi_);
    std::vector<double> y_airfoil(nxi_);
    
    // O-grid wraps around airfoil, starting from trailing edge
    // Go around ONCE (no wake cut!)
    for (int i = 0; i < nxi_; i++) {
        double theta = 2.0 * M_PI * i / (nxi_ - 1);  // Full circle: 0 → 2π
        
        // Cosine distribution for better clustering at TE/LE
        double s = 0.5 * (1.0 - std::cos(theta));  // 0 → 1 → 0 (goes around airfoil)
        
        // NACA 4-digit equations (same as C-grid)
        double yt = 5.0 * thick * (
            0.2969*std::sqrt(s) - 
            0.1260*s - 
            0.3516*s*s + 
            0.2843*s*s*s - 
            0.1015*s*s*s*s
        );
        
        double yc = 0.0;
        if (s < p) {
            yc = m / (p*p) * (2.0*p*s - s*s);
        } else {
            yc = m / ((1.0-p)*(1.0-p)) * ((1.0-2.0*p) + 2.0*p*s - s*s);
        }
        
        // Upper/lower surface based on which side of circle we're on
        double sign = (theta < M_PI) ? 1.0 : -1.0;  // Upper: 0→π, Lower: π→2π
        
        x_airfoil[i] = s * chord;
        y_airfoil[i] = (yc + sign * yt) * chord;
    }
    
    // STEP 2: Algebraic O-grid generation
    std::cout << "Step 2/3: Algebraic O-grid initialization...\n";
    
    // Use hyperbolic tangent stretching for better clustering near airfoil
    auto stretch_fcn = [](double eta, double beta) {
        // beta controls clustering: larger = more points near surface
        return std::tanh(beta * eta) / std::tanh(beta);
    };
    
    double beta = 4.0;  // Stretching parameter (higher = more clustering)
    
    for (int i = 0; i < nxi_; i++) {
        double theta = 2.0 * M_PI * i / (nxi_ - 1);
        
        // Airfoil point (inner boundary)
        double x_inner = x_airfoil[i];
        double y_inner = y_airfoil[i];
        
        // Find radial direction from airfoil point
        // Point from airfoil center toward this surface point
        double x_center = 0.5 * chord;
        double y_center = 0.0;
        
        double dx = x_inner - x_center;
        double dy = y_inner - y_center;
        double r_inner = std::sqrt(dx*dx + dy*dy);
        
        // Unit radial direction
        double nx = dx / r_inner;
        double ny = dy / r_inner;
        
        // March out in radial direction with stretching
        for (int j = 0; j < neta_; j++) {
            double eta_uniform = static_cast<double>(j) / (neta_ - 1);
            double eta = stretch_fcn(eta_uniform, beta);  // Clustered spacing
            
            // Radial distance from airfoil center
            double r = r_inner + eta * (farfield_radius - r_inner);
            
            // Position along radial line
            double x = x_center + r * nx;
            double y = y_center + r * ny;
            
            // Rotate by angle of attack
            double x_rot = std::cos(alpha) * x - std::sin(alpha) * y;
            double y_rot = std::sin(alpha) * x + std::cos(alpha) * y;
            
            int idx = j * nxi_ + i;
            x_[idx] = x_rot;
            y_[idx] = y_rot;
        }
    }
    
    std::cout << "  O-grid topology initialized (radial lines, no wake cut!)\n";
    
    // STEP 3: Elliptic smoothing (same as C-grid)
    std::cout << "Step 3/3: Elliptic smoothing...\n";
    smoothElliptic(10000, 5e-8);
    
    // STEP 4: Compute metrics
    computeMetrics();
    
    std::cout << "\nO-grid generation complete!\n\n";
}
