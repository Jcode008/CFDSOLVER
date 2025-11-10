#include "grid.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

Grid::Grid(int nx_, int ny_, double Lx_, double Ly_)
    : nx(nx_), ny(ny_), Lx(Lx_), Ly(Ly_)
{
    dx = Lx / nx;
    dy = Ly / ny;

    // Initialize mask to all fluid (1)
    mask = std::vector<std::vector<int>>(ny, std::vector<int>(nx, 1));
}

bool Grid::isFluid(int i, int j) const {
    if (i < 0 || i >= ny || j < 0 || j >= nx) {
        throw std::out_of_range("Grid::isFluid index out of range");
    }
    return mask[i][j] == 1;
}

// NACA 4-digit airfoil mask with rotation
void Grid::setAirfoilMask(double m, double p, double t, double x_pos, double y_pos, double chord, double angle)
{
    double angle_rad = angle * M_PI / 180.0;
    double cos_a = std::cos(angle_rad);
    double sin_a = std::sin(angle_rad);

    for(int i = 0; i < ny; i++){
        double y = i * dy;
        for(int j = 0; j < nx; j++){
            double x = j * dx;
            
            // Translate and rotate point to airfoil coordinates
            double dx_rel = x - x_pos;
            double dy_rel = y - y_pos;
            double x_local = dx_rel * cos_a + dy_rel * sin_a;
            double y_local = -dx_rel * sin_a + dy_rel * cos_a;
            
            // Check if point is within chord bounds
            if (x_local >= 0 && x_local <= chord) {
                double x_c = x_local / chord;
                
                // NACA 4-digit thickness distribution
                double yt = 5.0 * t * chord * (0.2969*std::sqrt(x_c) - 0.1260*x_c - 
                            0.3516*x_c*x_c + 0.2843*x_c*x_c*x_c - 0.1015*x_c*x_c*x_c*x_c);
                
                // Camber line (simplified)
                double yc = 0.0;
                if (p > 0) {
                    if (x_c < p) {
                        yc = m * chord * (2*p*x_c - x_c*x_c) / (p*p);
                    } else {
                        yc = m * chord * ((1-2*p) + 2*p*x_c - x_c*x_c) / ((1-p)*(1-p));
                    }
                }
                
                // Check if point is inside airfoil
                if (std::abs(y_local - yc) <= yt) {
                    mask[i][j] = 0; // mark as solid
                }
            }
        }
    }
}