#include "mesh_generator.hpp"
#include <iostream>

int main() {
    std::cout << "=========================================\n";
    std::cout << "  Body-Fitted O-Grid Generator Test\n";
    std::cout << "  (No wake cut - better than C-grid!)\n";
    std::cout << "=========================================\n\n";
    
    // Create grid generator
    // Start with medium resolution to test O-grid
    int nxi = 120;   // ξ direction (wraps around airfoil)
    int neta = 60;   // η direction (surface → farfield) - more points radially
    
    MeshGenerator mesh(nxi, neta);
    
    // Generate O-grid around NACA 2412
    double chord = 0.6;             // 0.6 m chord (same as your current setup)
    double alpha_deg = 0.0;         // 0° angle of attack (symmetric flow)
    double farfield_radius = 1.8;   // 3× chord (tighter for better aspect ratio)
    
    mesh.generateOGrid(chord, alpha_deg, farfield_radius);
    
    // Export for visualization
    std::cout << "Exporting grid files...\n";
    mesh.exportGrid("grid.csv");
    mesh.exportMetrics("metrics.csv");
    
    std::cout << "\n✓ Grid generation test complete!\n\n";
    std::cout << "Next steps:\n";
    std::cout << "  1. Visualize grid: python visualize_grid.py\n";
    std::cout << "  2. Check metrics (Jacobian should be >0 everywhere)\n";
    std::cout << "  3. If good, integrate with CFD solver\n\n";
    
    return 0;
}
