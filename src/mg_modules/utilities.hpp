#ifndef TYPES_HPP_
#define TYPES_HPP_

#include <vector>
#include <tuple>
#include <iostream>
#include <array>

namespace multigrid {

    enum CycleType{vCycle, wCycle, fCycle };
    enum BcType{dirichlet, neumann};
    struct BoundaryPoint {
        BcType bcType;
        double value;
    };
    typedef std::array<BoundaryPoint,1> Boundary;




};


#endif