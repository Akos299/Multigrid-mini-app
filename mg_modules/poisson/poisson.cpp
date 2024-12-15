

#include "poisson.hpp"

template<typename T, int Ndim>
Poisson<T,Ndim>::Poisson(const multigrid::Settings<T>& setting, std::array<double, 3> &starts, std::array<double, 3> &ends, multigrid::TransfertOperator &rest_operator,
                                          multigrid::TransfertOperator &inter_operator){

                                          }
template<typename T, int Ndim>
Poisson<T,Ndim>::~Poisson(){}

template<typename T, int Ndim>
void Poisson<T,Ndim>::solve()
{
    fmg(this->cycle_type);
}
        
