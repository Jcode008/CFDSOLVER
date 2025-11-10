#include "field.hpp"
#include <vector>

Field::Field(const Grid& grid, double U_infty) {
    nx_ = grid.nx;
    ny_ = grid.ny;
    size_t total = static_cast<size_t>(nx_) * ny_;
    u_.assign(total, U_infty);
    v_.assign(total, 0.0);
    p_.assign(total, 0.0);
    phi_.assign(total, 0.0);  // Pressure correction field
}

void Field::applyBoundaryConditions(const Grid& grid, double U_infty) {
    // Inlet
    for (int i = 0; i < grid.ny; ++i) {
        u(i, 0) = U_infty;
        v(i, 0) = 0.0;
    }

    // Outlet (zero-gradient)
    for (int i = 0; i < grid.ny; ++i) {
        u(i, grid.nx-1) = u(i, grid.nx-2);
        v(i, grid.nx-1) = v(i, grid.nx-2);
        p(i, grid.nx-1) = 0.0;
    }

    // Top/bottom (free-slip)
    for (int j = 0; j < grid.nx; ++j) {
        u(0, j) = u(1, j);
        u(grid.ny-1, j) = u(grid.ny-2, j);
        v(0, j) = 0.0;
        v(grid.ny-1, j) = 0.0;
    }
    // Airfoil solid cells
    for (int i = 0; i < grid.ny; ++i) {
        for (int j = 0; j < grid.nx; ++j) {
            if (!grid.isFluid(i,j)) {
                u(i, j) = 0.0;
                v(i, j) = 0.0;
            }
        }
    }
}
