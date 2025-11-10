#include "curvilinear_field.hpp"

CurvilinearField::CurvilinearField(int nxi, int neta) 
    : nxi_(nxi), neta_(neta) {
    int total = nxi_ * neta_;
    u_.resize(total, 0.0);
    v_.resize(total, 0.0);
    p_.resize(total, 0.0);
    u_star_.resize(total, 0.0);
    v_star_.resize(total, 0.0);
}
