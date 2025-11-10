#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <vector>
#include "field.hpp"
#include "grid.hpp"

//solver class handles the Navier-Stokes Computation
class Solver {
public:
    Solver(Grid &grid, Field &field, double rho, double nu, double dt);

    //Core Steps
    void computeIntermediateVelocities();
    void solvePressurePoisson(int nit = 300, double tol = 1e-4, double omega = 1.7);
    void correctVelocities();

    //Time Stepping 
    void step(); //single time step 

private:
    Grid &grid;
    Field &field;
    double rho; //density
    double nu;  //kinematic viscosity
    double dt;  //time step size

    //Temporary storage
    std::vector<double> u_star;
    std::vector<double> v_star;
    size_t totalCells = 0;

    inline size_t idx(int i, int j) const { return static_cast<size_t>(i) * grid.nx + j; }
};

#endif