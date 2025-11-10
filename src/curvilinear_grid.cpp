#include "curvilinear_grid.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>

CurvilinearGrid::CurvilinearGrid(int nxi, int neta) 
    : nxi_(nxi), neta_(neta) {
    int total = nxi_ * neta_;
    x_.resize(total);
    y_.resize(total);
    J_.resize(total);
    xi_x_.resize(total);
    xi_y_.resize(total);
    eta_x_.resize(total);
    eta_y_.resize(total);
}

void CurvilinearGrid::loadFromFile(const std::string& grid_file, const std::string& metrics_file) {
    std::cout << "Loading curvilinear grid from files...\n";
    
    // Load grid coordinates
    std::ifstream grid_in(grid_file);
    if (!grid_in) {
        throw std::runtime_error("Cannot open grid file: " + grid_file);
    }
    
    std::string header;
    std::getline(grid_in, header);  // Skip header: x,y,i,j
    
    int points_loaded = 0;
    std::string line;
    while (std::getline(grid_in, line)) {
        std::istringstream iss(line);
        double x, y;
        int i, j;
        char comma;
        
        if (iss >> x >> comma >> y >> comma >> i >> comma >> j) {
            int idx = j * nxi_ + i;
            x_[idx] = x;
            y_[idx] = y;
            points_loaded++;
        }
    }
    grid_in.close();
    
    std::cout << "  Loaded " << points_loaded << " grid points\n";
    
    // Load metrics
    std::ifstream metrics_in(metrics_file);
    if (!metrics_in) {
        throw std::runtime_error("Cannot open metrics file: " + metrics_file);
    }
    
    std::getline(metrics_in, header);  // Skip header
    
    int metrics_loaded = 0;
    while (std::getline(metrics_in, line)) {
        std::istringstream iss(line);
        int i, j;
        double J, xi_x, xi_y, eta_x, eta_y;
        char comma;
        
        if (iss >> i >> comma >> j >> comma >> J >> comma 
                >> xi_x >> comma >> xi_y >> comma 
                >> eta_x >> comma >> eta_y) {
            int idx = j * nxi_ + i;
            J_[idx] = J;
            xi_x_[idx] = xi_x;
            xi_y_[idx] = xi_y;
            eta_x_[idx] = eta_x;
            eta_y_[idx] = eta_y;
            metrics_loaded++;
        }
    }
    metrics_in.close();
    
    std::cout << "  Loaded " << metrics_loaded << " metric values\n";
    
    // Validate grid
    int negative_jacobians = 0;
    int zero_jacobians = 0;
    double min_J = 1e10, max_J = -1e10;
    
    for (int j = 0; j < neta_; j++) {
        for (int i = 0; i < nxi_; i++) {
            double J_val = J(i, j);
            if (J_val < 0.0) negative_jacobians++;
            if (J_val == 0.0) zero_jacobians++;
            if (J_val > 0.0) {  // Only track positive values
                min_J = std::min(min_J, J_val);
                max_J = std::max(max_J, J_val);
            }
        }
    }
    
    std::cout << "  Jacobian range: [" << min_J << ", " << max_J << "]\n";
    
    if (negative_jacobians > 0) {
        std::cout << "  ⚠ WARNING: " << negative_jacobians 
                  << " points have negative Jacobian!\n";
    }
    if (zero_jacobians > 0) {
        std::cout << "  ⚠ WARNING: " << zero_jacobians 
                  << " points have zero Jacobian (boundary points not exported by grid generator)\n";
    }
    if (negative_jacobians == 0 && zero_jacobians == 0) {
        std::cout << "  ✓ All Jacobians positive - valid grid!\n";
    }
    
    std::cout << "\n";
}

double CurvilinearGrid::cellArea(int i, int j) const {
    // Cell area in physical space = J * Δξ * Δη
    // For uniform computational grid: Δξ = Δη = 1/(n-1)
    double dxi = 1.0 / (nxi_ - 1);
    double deta = 1.0 / (neta_ - 1);
    return J(i, j) * dxi * deta;
}
