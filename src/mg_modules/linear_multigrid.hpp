#ifndef __LINEAR_MULTIGRID__
#define __LINEAR_MULTIGRID__

#include "multigrid_base.hpp"

namespace multigrid {

template<typename  T, int Ndim>
class LinearMultigrid : public MultigridBase<T, Ndim>
{
    LinearMultigrid(const Settings<T>& settings): MultigridBase<T, Ndim>::MultigridBase(settings)
    {};
    virtual ~LinearMultigrid(){};

    virtual void multigrid();

};

}


#endif 