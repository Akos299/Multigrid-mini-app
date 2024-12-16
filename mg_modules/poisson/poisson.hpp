#ifndef __POISSON_HPP__
#define __POISSON_HPP__

#include "../include/multigrid_base.hpp"

template<typename T, int Ndim>
class Poisson: public multigrid::MultigridBase<T,Ndim>{
public:
    Poisson(const multigrid::Settings<T>& settings, std::array<double, 3> &starts, std::array<double, 3> &ends, multigrid::TransfertOperator &rest_operator,
                                          multigrid::TransfertOperator &inter_operator);
    virtual ~Poisson();
    void solve();

};

#endif

