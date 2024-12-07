#ifndef TYPES_HPP_
#define TYPES_HPP_

#include <vector>
#include <tuple>
#include <iostream>
#include <array>
#include <math.h>

namespace multigrid {

    enum CycleType{vCycle, wCycle, fCycle,userdef };
    enum BcType{dirichlet, neumann};
    enum TransfertOperator{PointAverage, Linear, Khalil, Kwak};
    struct BoundaryPoint {
        BcType bcType;
        double value;
    };
    typedef std::array<BoundaryPoint,1> Boundary;

/** @brief // This code is contributed by Aditya Kumar (adityakumar129). From Geeksforgeeks
 */
int power(int x, int n)
{
    // If x^0 return 1
    if (n == 0)
        return 1;
    // If we need to find of 0^y
    if (x == 0)
        return 0;
    // For all other cases
    return x * power(x, n - 1);
}

typedef struct {
    double dx, dy, dz, dxx, dyy, dzz, dxy, dxz, dyz;
} Derivatives;


};


#endif