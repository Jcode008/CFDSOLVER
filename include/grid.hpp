#ifndef GRID_HPP
#define GRID_HPP

#include <vector>
#include <cmath>

class Grid { 
    public: 
    int nx, ny;   //number of cells in x and y 
    double Lx, Ly; //physical Domain dimensions
    double dx, dy; //grid spacing
    std::vector<std::vector<int>> mask; //1:fluid, 0:solid

    Grid(int nx_, int ny_, double Lx_, double Ly_);

    bool isFluid(int i, int j)const; 

    void setAirfoilMask(double m, double p, double t, double x_pos, double y_pos, double chord, double angle);
};

#endif

