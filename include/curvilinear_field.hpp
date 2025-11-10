#pragma once
#include <vector>

// Field variables on curvilinear grid
// Instead of u(x,y), we store u(ξ,η)
// The physical meaning is the same, just indexed differently

class CurvilinearField {
public:
    CurvilinearField(int nxi, int neta);
    
    // Access field values at computational grid point (i,j)
    double& u(int i, int j) { return u_[j * nxi_ + i]; }
    double& v(int i, int j) { return v_[j * nxi_ + i]; }
    double& p(int i, int j) { return p_[j * nxi_ + i]; }
    
    const double& u(int i, int j) const { return u_[j * nxi_ + i]; }
    const double& v(int i, int j) const { return v_[j * nxi_ + i]; }
    const double& p(int i, int j) const { return p_[j * nxi_ + i]; }
    
    // Intermediate velocities (for fractional step)
    double& u_star(int i, int j) { return u_star_[j * nxi_ + i]; }
    double& v_star(int i, int j) { return v_star_[j * nxi_ + i]; }
    
    const double& u_star(int i, int j) const { return u_star_[j * nxi_ + i]; }
    const double& v_star(int i, int j) const { return v_star_[j * nxi_ + i]; }
    
    int nxi() const { return nxi_; }
    int neta() const { return neta_; }
    
private:
    int nxi_, neta_;
    std::vector<double> u_, v_, p_;
    std::vector<double> u_star_, v_star_;
};
