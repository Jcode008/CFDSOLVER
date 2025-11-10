#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include "grid.hpp"
#include "field.hpp"
#include "solver.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
/**
 * @brief Main function for the CFD solver simulation.
 * 
 * This function sets up and initializes a CFD simulation with the following steps:
 * 
 * 1. Domain and Grid Setup:
 *    - Lx: Domain length in x-direction (meters)
 *    - Ly: Domain length in y-direction (meters)
 *    - nx: Number of grid points in x-direction
 *    - ny: Number of grid points in y-direction
 * 
 * 2. Fluid Properties:
 *    - rho: Fluid density (kg/m³)
 *    - nu: Kinematic viscosity (m²/s)
 *    - U_infty: Freestream velocity (m/s)
 * 
 * 3. Field Initialization:
 *    - Initializes velocity field with u = U_infty, v = 0, and pressure p = 0
 * 
 * 4. Airfoil Mask (NACA 2412):
 *    - Parameter 1 (0.02): Maximum camber as fraction of chord
 *    - Parameter 2 (0.4): Position of maximum camber as fraction of chord
 *    - Parameter 3 (0.12): Maximum thickness as fraction of chord
 *    - Parameter 4 (0.25*Lx): Airfoil x-position (center/leading edge)
 *    - Parameter 5 (0.5*Ly): Airfoil y-position (vertical center)
 *    - Parameter 6 (0.3*Lx): Airfoil chord length
 *    - Parameter 7 (-10.0): Angle of attack in degrees
 * 
 */

int main() {
    double Lx = 4.0, Ly = 2.0; //Domain size - smaller, taller domain
    int nx = 800, ny = 400;    // DOUBLED resolution for better airfoil resolution
    Grid grid(nx, ny, Lx, Ly);


    //-------------------------
    // 2. Fluid Properties
    //-------------------------
    double rho = 1.0;    // Density (kg/m³)
    double nu  = 2e-4;   // Kinematic viscosity - Re = U*c/nu = 5*0.4/2e-4 = 10,000
    double U_infty = 5.0; // Freestream velocity (m/s)

    //-------------------------
    // 3. Field Initialization
    //--------------------------
    Field field(grid, U_infty);


    //-------------------------
    // 4. NACA Airfoil Mask
    //--------------------------
    
    // NACA 2412: Classic cambered airfoil
    // 2 = 2% maximum camber
    // 4 = position of max camber at 40% chord
    // 12 = 12% maximum thickness
    double m = 0.02;      // Max camber (2%)
    double p = 0.4;       // Position of max camber (40% chord)
    double t = 0.12;      // Max thickness (12% chord)
    double chord = 0.6;   // 0.6m chord length (60 grid cells)
    double x_le = 0.6;    // Leading edge at x=0.6m
    double y_c = 1.0;     // Centered vertically at y=1.0m
    double alpha = 0.0;   // Angle of attack (degrees) - 5° for good lift
    
    // Create airfoil using existing NACA function
    grid.setAirfoilMask(m, p, t, x_le, y_c, chord, alpha);
    
    // Count solid cells for verification
    int solid_count = 0;
    for(int i = 0; i < ny; i++){
        for(int j = 0; j < nx; j++){
            if(!grid.isFluid(i, j)) solid_count++;
        }
    }
    
    std::cout << "NACA " << int(m*100) << int(p*10) << int(t*100) 
              << " airfoil created at alpha=" << alpha << "°" << std::endl;
    std::cout << "Chord = " << chord << "m, Reynolds number = " 
              << U_infty*chord/nu << std::endl;
    std::cout << "Solid cells: " << solid_count << " / " << (nx*ny) << std::endl;

    //-------------------------
    // 5. Solver Initialization
    //--------------------------
    double dt = 1.25e-5; // QUARTERED timestep for finer grid (CFL stability)
    Solver solver(grid, field, rho, nu, dt);

    //-------------------------
    // 6. Time Stepping Loop
    //-------------------------
    int nt = 2500; // Reduced timesteps for finer grid (same physical time)
    int snapshotInterval = 100; // Save every 100 steps (same number of snapshots)

    for(int n = 0; n<nt; n++){
        if (n % 50 == 0) std::cout << "Starting timestep " << n << std::endl;
        
        // inlet BC
        for(int i =0 ; i < ny; i++){
            field.u(i, 0) = U_infty;
            field.v(i, 0) = 0.0;
        }

        // top/bottom BC (free-slip)
        for(int j =0 ; j < nx; j++){
            field.u(0, j) = field.u(1, j);
            field.u(ny-1, j) = field.u(ny-2, j);
            field.v(0, j) = 0.0;
            field.v(ny-1, j) = 0.0;
        }

        // outlet BC: zero-gradient with sponge blending
        for (int i = 0; i < ny; ++i) {
            field.u(i, nx-1) = field.u(i, nx-2);
            field.v(i, nx-1) = field.v(i, nx-2);
            field.p(i, nx-1) = field.p(i, nx-2);
        }

        int spongeWidth = 80;
        int startCol = nx - spongeWidth;
        if (startCol < 1) startCol = 1;
        for (int j = startCol; j < nx; ++j) {
            double w = static_cast<double>(j - startCol + 1) / spongeWidth;
            for (int i = 0; i < ny; ++i) {
                if (!grid.isFluid(i, j)) continue;
                field.u(i, j) = (1.0 - w) * field.u(i, j) + w * U_infty;
                field.v(i, j) *= (1.0 - w);
            }
        }

        //solid Cells BC
         for(int i=0; i<ny; i++){
            for(int j=0; j<nx; j++){
                if(!grid.isFluid(i,j)){
                    field.u(i, j) = 0.0;
                    field.v(i, j) = 0.0;
                }
            }
        }
        
        //solver step 
        solver.step();
        
        if (n % 50 == 0) std::cout << "Completed solver step for " << n << std::endl;
        
        if(n % 50 == 0){
            double mass_in = 0.0;
            double mass_out = 0.0;
            for(int i = 0; i < ny; i++){
                mass_in += field.u(i, 0);
                mass_out += field.u(i, nx-1);
            }
            std::cout << "Step " << n << ": Mass in=" << mass_in << ", Mass out=" << mass_out << std::endl;
        }

        //Save Snapshots 
        if(n % snapshotInterval == 0){
            // Save u velocity
            std::string filename_u = "u_" + std::to_string(n) + ".csv";
            std::ofstream fout_u(filename_u);
            for(int i = 0; i < ny; i++){
                for(int j = 0; j < nx; j++){
                    if(grid.isFluid(i, j)){
                        fout_u << field.u(i, j);
                    } else {
                        fout_u << "nan";
                    }
                    if(j < nx -1) fout_u << ",";
                }
                fout_u << "\n";
            }
            fout_u.close();
            
            // Save v velocity
            std::string filename_v = "v_" + std::to_string(n) + ".csv";
            std::ofstream fout_v(filename_v);
            for(int i = 0; i < ny; i++){
                for(int j = 0; j < nx; j++){
                    if(grid.isFluid(i, j)){
                        fout_v << field.v(i, j);
                    } else {
                        fout_v << "nan";
                    }
                    if(j < nx -1) fout_v << ",";
                }
                fout_v << "\n";
            }
            fout_v.close();
            
            // Save pressure
            std::string filename_p = "p_" + std::to_string(n) + ".csv";
            std::ofstream fout_p(filename_p);
            for(int i = 0; i < ny; i++){
                for(int j = 0; j < nx; j++){
                    if(grid.isFluid(i, j)){
                        fout_p << field.p(i, j);
                    } else {
                        fout_p << "nan";
                    }
                    if(j < nx -1) fout_p << ",";
                }
                fout_p << "\n";
            }
            fout_p.close();
        }
        if(n % 50 == 0){
            std::cout << "Completed time step: " << n << "/" << nt << std::endl;
        }


    }

    std::cout << "Simulation complete." << std::endl;
    return 0;
}