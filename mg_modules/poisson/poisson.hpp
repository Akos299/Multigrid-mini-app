#ifndef __POISSON_HPP__
#define __POISSON_HPP__

#include "../include/multigrid_base.hpp"
#include "../include/utilities.hpp"
#include <array>
#include "../include/ndArray.h"


namespace multigrid{

// template<class T, int Ndim>
// class Poisson: public MultigridBase<T,Ndim>{
// public:
//     // Poisson(const Settings<T>& settings, std::array<double, 3> &starts, std::array<double, 3> &ends, TransfertOperator &rest_operator,
//     //                                       TransfertOperator &inter_operator): MultigridBase<T,Ndim>(settings, starts, ends, rest_operator, inter_operator)
//     //                                       {
//     //                                             std::cout << "constructor" << "\n\n";
//     //                                       };

//     Poisson(my_settings& settings): MultigridBase<T, Ndim>::MultigridBase(settings){
//         /***/
//     };
//     virtual ~Poisson(){};
// //    virtual void solve();

// };

template<class T, int Ndim>
class Poisson : public MultigridBase<T,Ndim> {
public:
    // Poisson(const Settings<T>& settings, std::array<double, 3> &starts, std::array<double, 3> &ends, TransfertOperator &rest_operator,
    //                                       TransfertOperator &inter_operator): MultigridBase<T,Ndim>(settings, starts, ends, rest_operator, inter_operator)
    //                                       {
    //                                             std::cout << "constructor" << "\n\n";
    //                                       };

    Poisson(my_settings& settings,std::array<double, 3> &starts, std::array<double, 3> &ends, TransfertOperator &rest_operator,
                                          TransfertOperator &inter_operator) : MultigridBase<T,Ndim>::MultigridBase(settings,starts, ends, rest_operator,inter_operator){
        /***/ std::cout << "I'm in the poisson constructor" << "\n";
    };
    virtual ~Poisson(){};
   virtual void solve();

};

template<typename T, int Ndim>
void multigrid::Poisson<T,Ndim>::solve()
{
    this->fmg(this->cycle_type);
}

}

#endif

