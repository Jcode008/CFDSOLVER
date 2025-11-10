#ifndef FIELD_HPP
#define FIELD_HPP

#include "grid.hpp"
#include <vector>

class Field {
public:
    Field(const Grid& grid, double U_infty);

    void applyBoundaryConditions(const Grid& grid, double U_infty);

    inline double& u(int i, int j) { return u_[index(i, j)]; }
    inline double& v(int i, int j) { return v_[index(i, j)]; }
    inline double& p(int i, int j) { return p_[index(i, j)]; }
    inline double& phi(int i, int j) { return phi_[index(i, j)]; }

    inline const double& u(int i, int j) const { return u_[index(i, j)]; }
    inline const double& v(int i, int j) const { return v_[index(i, j)]; }
    inline const double& p(int i, int j) const { return p_[index(i, j)]; }
    inline const double& phi(int i, int j) const { return phi_[index(i, j)]; }

    inline int nx() const { return nx_; }
    inline int ny() const { return ny_; }

    inline std::vector<double>& rawU() { return u_; }
    inline std::vector<double>& rawV() { return v_; }
    inline std::vector<double>& rawP() { return p_; }
    inline std::vector<double>& rawPhi() { return phi_; }
    inline const std::vector<double>& rawU() const { return u_; }
    inline const std::vector<double>& rawV() const { return v_; }
    inline const std::vector<double>& rawP() const { return p_; }
    inline const std::vector<double>& rawPhi() const { return phi_; }

private:
    int nx_ = 0;
    int ny_ = 0;
    std::vector<double> u_;
    std::vector<double> v_;
    std::vector<double> p_;
    std::vector<double> phi_;  // Pressure correction

    inline size_t index(int i, int j) const { return static_cast<size_t>(i) * nx_ + j; }
};

#endif
