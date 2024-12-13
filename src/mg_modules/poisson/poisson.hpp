#ifndef __POISSON__
#define __POISSON__

#include "../multigrid_base.hpp"

template<typename, int Ndim>
class Poisson: public MultigridBase<T,Ndim>{
public:
    Poisson(const multigrid::Settings<T>& settings, std::array<double, 3> &starts, std::array<double, 3> &ends, multigrid::TransfertOperator &rest_operator,
                                          multigrid::TransfertOperator &inter_operator);
    virtual ~Poisson();
    void solve();


};


#endif