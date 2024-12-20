#ifndef __POISSON_HPP__
#define __POISSON_HPP__

#include "multigrid_base.hpp"
#include "ndArray.h"
#include "utilities.hpp"
#include <array>
#include <cstddef>
#include <string>

namespace multigrid {

// template<class T, int Ndim>
// class Poisson: public MultigridBase<T,Ndim>{
// public:
//     // Poisson(const Settings<T>& settings, std::array<double, 3> &starts,
//     std::array<double, 3> &ends, TransfertOperator &rest_operator,
//     //                                       TransfertOperator
//     &inter_operator): MultigridBase<T,Ndim>(settings, starts, ends,
//     rest_operator, inter_operator)
//     //                                       {
//     //                                             std::cout << "constructor"
//     << "\n\n";
//     //                                       };

//     Poisson(my_settings& settings): MultigridBase<T,
//     Ndim>::MultigridBase(settings){
//         /***/
//     };
//     virtual ~Poisson(){};
// //    virtual void solve();

// };

template <class T, int Ndim> class Poisson : public MultigridBase<T, Ndim> {
public:
  // Poisson(const Settings<T>& settings, std::array<double, 3> &starts,
  // std::array<double, 3> &ends, TransfertOperator &rest_operator,
  //                                       TransfertOperator &inter_operator):
  //                                       MultigridBase<T,Ndim>(settings,
  //                                       starts, ends, rest_operator,
  //                                       inter_operator)
  //                                       {
  //                                             std::cout << "constructor" <<
  //                                             "\n\n";
  //                                       };

  Poisson(my_settings &settings, std::array<double, 3> &starts,
          std::array<double, 3> &ends, TransfertOperator &rest_operator,
          TransfertOperator &inter_operator)
      : MultigridBase<T, Ndim>::MultigridBase(settings, starts, ends,
                                              rest_operator, inter_operator) {
    /***/ std::cout << "I'm in the poisson constructor" << "\n";
  };
  virtual ~Poisson() {};
  virtual void solve();
  virtual void set_level_solution(size_t level, T value) ;
  virtual void set_level_source(size_t level, T value);
  virtual void print_level_data(size_t level, std::string data_value);
};

template <typename T, int Ndim> void multigrid::Poisson<T, Ndim>::solve() {
  this->fmg(this->cycle_type);
}

/**/
template<typename T, int Ndim>
void Poisson<T, Ndim>::set_level_source(size_t level, T value) 
{
    std::cout << "Poisson::set_level_source::start \n";
    MultigridBase<T,Ndim>::set_level_source(level, value);
    std::cout << "Poisson::set_level_source::end \n";
}

/**/
template<typename T, int Ndim>
void Poisson<T, Ndim>::set_level_solution(size_t level, T value)
{
    std::cout << "Poisson::set_level_solution::start \n";
    MultigridBase<T,Ndim>::set_level_solution(level, value);
    std::cout << "Poisson::set_level_solution::end \n";
}

/**/
template<typename T, int Ndim>
void Poisson<T, Ndim>::print_level_data(size_t level, std::string data_value)
{
    std::cout << "Poisson::print_level_data::start \n";
    MultigridBase<T, Ndim>::print_level_data(level, data_value);
    std::cout << "Poisson::print_level_data::end \n";
}



 



} // namespace multigrid

#endif
