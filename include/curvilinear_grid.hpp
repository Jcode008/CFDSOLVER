#pragma once
#include <vector>
#include <string>

// Curvilinear body-fitted grid for CFD solver
// Stores physical coordinates x(ξ,η), y(ξ,η) and transformation metrics
// 
// WHY THIS IS DIFFERENT FROM CARTESIAN:
// - Grid points don't align with x,y axes
// - Grid spacing non-uniform in physical space
// - Need metrics (ξ_x, ξ_y, η_x, η_y, J) to transform derivatives
//
// Reference: Thompson et al. (1985) Ch. 5

class CurvilinearGrid {
public:
    CurvilinearGrid(int nxi, int neta);
    
    // Load grid from mesh generator output
    void loadFromFile(const std::string& grid_file, const std::string& metrics_file);
    
    // Grid dimensions in computational space
    int nxi() const { return nxi_; }
    int neta() const { return neta_; }
    
    // Physical coordinates at grid point (i,j)
    double x(int i, int j) const { return x_[j * nxi_ + i]; }
    double y(int i, int j) const { return y_[j * nxi_ + i]; }
    
    // Transformation metrics (Thompson Eq. 5.11)
    // These are needed to transform ∂/∂x, ∂/∂y to ∂/∂ξ, ∂/∂η
    double J(int i, int j) const { return J_[j * nxi_ + i]; }       // Jacobian
    double xi_x(int i, int j) const { return xi_x_[j * nxi_ + i]; } // ∂ξ/∂x
    double xi_y(int i, int j) const { return xi_y_[j * nxi_ + i]; } // ∂ξ/∂y
    double eta_x(int i, int j) const { return eta_x_[j * nxi_ + i]; } // ∂η/∂x
    double eta_y(int i, int j) const { return eta_y_[j * nxi_ + i]; } // ∂η/∂y
    
    // Cell volume in physical space (for divergence calculations)
    // In 2D: area = J * Δξ * Δη
    double cellArea(int i, int j) const;
    
    // Check if point (i,j) is on a boundary
    bool isAirfoil(int i, int j) const { return j == 0; }
    bool isFarfield(int i, int j) const { return j == neta_ - 1; }
    bool isWakeCut(int i, int j) const { return i == 0 || i == nxi_ - 1; }
    
private:
    int nxi_;   // Grid points in ξ direction (around airfoil)
    int neta_;  // Grid points in η direction (normal to airfoil)
    
    // Physical coordinates
    std::vector<double> x_, y_;
    
    // Transformation metrics
    std::vector<double> J_;      // Jacobian = x_ξ*y_η - x_η*y_ξ
    std::vector<double> xi_x_;   // ∂ξ/∂x = y_η/J
    std::vector<double> xi_y_;   // ∂ξ/∂y = -x_η/J
    std::vector<double> eta_x_;  // ∂η/∂x = -y_ξ/J
    std::vector<double> eta_y_;  // ∂η/∂y = x_ξ/J
};
