#include "curvilinear_grid.hpp"
#include "curvilinear_field.hpp"
#include "curvilinear_solver.hpp"
#include "io.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Helper function to export field data
void exportFieldCSV(const std::string& filename, const CurvilinearField& field, 
                    int nxi, int neta, char var_type) {
    std::ofstream file(filename);
    file << "i,j,value\n";
    for (int j = 0; j < neta; j++) {
        for (int i = 0; i < nxi; i++) {
            double value = 0.0;
            if (var_type == 'u') value = field.u(i, j);
            else if (var_type == 'v') value = field.v(i, j);
            else if (var_type == 'p') value = field.p(i, j);
            file << i << "," << j << "," << value << "\n";
        }
    }
    file.close();
}

int main() {
    std::cout << "========================================\n";
    std::cout << "  Curvilinear CFD Solver Test\n";
    std::cout << "  NACA 2412 at α=0° (Zero AoA)\n";
    std::cout << "========================================\n\n";
    
    // STEP 1: Load body-fitted grid
    std::cout << "STEP 1: Loading grid...\n";
    CurvilinearGrid grid(120, 60);  // O-grid resolution (updated from 40 to 60)
    
    try {
        grid.loadFromFile("grid.csv", "metrics.csv");
    } catch (const std::exception& e) {
        std::cerr << "Error loading grid: " << e.what() << "\n";
        std::cerr << "\nMake sure to run TestMesh first to generate grid!\n";
        std::cerr << "  cd build\n";
        std::cerr << "  .\\Release\\TestMesh.exe\n\n";
        return 1;
    }
    
    // STEP 2: Initialize flow field
    std::cout << "STEP 2: Initializing flow field...\n";
    CurvilinearField field(grid.nxi(), grid.neta());
    
    // Set initial conditions: uniform horizontal freestream
    // IMPORTANT: Grid has airfoil at α=-5° (nose down rotation)
    // With horizontal flow, this gives POSITIVE α=+5° angle of attack!
    // (Flow hits airfoil from below when airfoil is rotated nose-down)
    double U_inf = 5.0;  // m/s horizontal
    
    for (int j = 0; j < grid.neta(); j++) {
        for (int i = 0; i < grid.nxi(); i++) {
            field.u(i, j) = U_inf;  // Horizontal flow
            field.v(i, j) = 0.0;     // No vertical component
            field.p(i, j) = 0.0;
        }
    }
    std::cout << "  Freestream: U∞ = " << U_inf << " m/s (horizontal)\n";
    std::cout << "  Grid rotation: α_grid = -5° → Effective α = +5°\n\n";
    
    // STEP 3: Setup solver
    std::cout << "STEP 3: Setting up solver...\n";
    double chord = 0.6;  // m
    double Re = 15000.0;
    double dt = 1e-5;    // Start conservative with coarser grid
    
    CurvilinearSolver solver(grid, Re, dt);
    
    // STEP 4: Run simulation
    std::cout << "STEP 4: Running simulation...\n";
    int n_steps = 2000;    // Run much longer to test stability
    int print_interval = 100;
    
    std::cout << "\nTimestep,MaxU,MinU,MaxV,MinV,MaxP,MinP\n";
    
    for (int step = 0; step <= n_steps; step++) {
        // Monitor solution
        if (step % print_interval == 0) {
            double max_u = -1e10, min_u = 1e10;
            double max_v = -1e10, min_v = 1e10;
            double max_p = -1e10, min_p = 1e10;
            bool has_nan = false;
            
            for (int j = 0; j < grid.neta(); j++) {
                for (int i = 0; i < grid.nxi(); i++) {
                    double u_val = field.u(i,j);
                    double v_val = field.v(i,j);
                    double p_val = field.p(i,j);
                    
                    if (std::isnan(u_val) || std::isnan(v_val) || std::isnan(p_val)) {
                        has_nan = true;
                        if (step == print_interval) {  // Print first occurrence
                            std::cout << "\n⚠ NaN at step " << step << ", i=" << i << ", j=" << j << "\n";
                            std::cout << "  u=" << u_val << ", v=" << v_val << ", p=" << p_val << "\n";
                        }
                    }
                    
                    max_u = std::max(max_u, u_val);
                    min_u = std::min(min_u, u_val);
                    max_v = std::max(max_v, v_val);
                    min_v = std::min(min_v, v_val);
                    max_p = std::max(max_p, p_val);
                    min_p = std::min(min_p, p_val);
                }
            }
            
            std::cout << step << ","
                      << max_u << "," << min_u << ","
                      << max_v << "," << min_v << ","
                      << max_p << "," << min_p << "\n";
            
            // Check for NaN
            if (has_nan) {
                std::cout << "\n⚠ ERROR: NaN detected! Simulation unstable.\n";
                std::cout << "Grid aspect ratio too high for explicit advection.\n";
                std::cout << "Need smaller timestep or fully implicit solver.\n";
                return 1;
            }
        }
        
        // Take timestep (AFTER printing, so step N prints state N, then advances to N+1)
        solver.step(field);
    }
    
    std::cout << "\n✓ Simulation completed successfully!\n";
    
    // Export final solution for visualization
    std::cout << "\nExporting solution data...\n";
    exportFieldCSV("u_final.csv", field, grid.nxi(), grid.neta(), 'u');
    exportFieldCSV("v_final.csv", field, grid.nxi(), grid.neta(), 'v');
    exportFieldCSV("p_final.csv", field, grid.nxi(), grid.neta(), 'p');
    std::cout << "  ✓ Exported u_final.csv, v_final.csv, p_final.csv\n";
    
    std::cout << "\nNext steps:\n";
    std::cout << "  1. Visualize flow: python analysis/visualize_flow.py\n";
    std::cout << "  2. Compute aerodynamic forces (lift, drag)\n";
    std::cout << "  3. Run parametric sweep (different α)\n\n";
    
    return 0;
}
