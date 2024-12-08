#ifndef TYPES_HPP_
#define TYPES_HPP_

#include <vector>
#include <tuple>
#include <iostream>
#include <array>
#include <math.h>
#include "ndArray.h"

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

typedef struct {
std::array<std::array<int, 3>, 3> Laplace2D  {{{0,1,0},{1,-6,1},{0,1,0}}};
} DiscreteOperator2D_Poisson;

typedef struct {
    std::array<std::array<int, 3>, 3> Laplace3D_L  {{{0,0,0},{0,1,0},{0,0,0}}}; // k-1
    std::array<std::array<int, 3>, 3> Laplace3D_C  {{{0,1,0},{1,-6,1},{0,1,0}}}; // k
    std::array<std::array<int, 3>, 3> Laplace3D_R  {{{0,0,0},{0,1,0},{0,0,0}}}; // k+1
} DiscreteOperator3D_Poisson;

enum DiscreteOperator{
    laplacian
};




/********************** Setting structure *****************
*   This is a structure to hold initial configuration of the multigrid solver. That's all the different components of the multigrid solver.
*/

    /** @brief structure to set initial parameter for the multigrid class */
    template <typename T>
    struct Settings
    {
        // T aspectRatio{1};
        size_t numberOfGrids{8}; // maxLevel
        size_t x_minimumResolution{4};
        size_t y_minimumResolution{4};
        size_t z_minimumResolution{4};
        T residualTolerance{1e-10};
        size_t maxIterations{400};
        CycleType cycleType{multigrid::wCycle};
        size_t Npre{1}, Npost{1}, maxLevel{21}, Nvcyle_cvg{20};
        T w{1.0};
    };


}


#endif