#include "../include/mgrid_stack.hpp"

namespace multigrid {

    template <typename T,int Ndim>
    template<TransfertOperator Op>
    inline void Stack<T, Ndim>::coarsen(size_t level)
    {
        restriction_operator<T, Op,Ndim>((*this)[level - 1], (*this)[level]);
    }

    template <typename T, int Ndim>
     template<TransfertOperator Op>
    inline void Stack<T, Ndim>::coarsen(size_t level, nd::ndArray<T, Ndim> &coarse_level_data)
    {
        restriction_operator<T,Op,Ndim>(coarse_level_data, (*this)[level]);
    }

    template <typename T, int Ndim>
     template<TransfertOperator Op>
    inline void Stack<T, Ndim>::refine(size_t level)
    {
        prolongation_operator<T,Op,Ndim>((*this)[level], (*this)[level + 1]);
    }

    template <typename T, int Ndim>
    template<TransfertOperator Op>
    inline void Stack<T, Ndim>::refine(size_t level, nd::ndArray<T, Ndim> &fine_level_data)
    {
        prolongation_operator<T,Op,Ndim>((*this)[level], fine_level_data);
    }

}