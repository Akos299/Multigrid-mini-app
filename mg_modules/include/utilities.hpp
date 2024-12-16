#ifndef __UTILITIES_HPP__
#define __UTILITIES_HPP__

#include <array>
#include <math.h>

namespace multigrid {

    /** @brief Different type of multigrid iteration cycle */
    enum CycleType{vCycle, wCycle, fCycle,userdef };

    /** @brief Different type of Boundary condition*/
    enum BcType{Periodic, ZeroGradient, ZeroFixed, Isolated, UserDef};

    /** @brief Different type of Transfert Operator*/
    enum TransfertOperator{PointAverage, Linear, Khalil, Kwak,TreeCubic};

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


/** @brief a structure for derivatives */
typedef struct {
    double dx, dy, dz, dxx, dyy, dzz, dxy, dxz, dyz;
} Derivatives;


// /** @brief a structure for Laplacian in 2D*/
// typedef struct {
// std::array<std::array<int, 3>, 3> Laplace2D  = {{{0,1,0},{1,-6,1},{0,1,0}}};
// } DiscreteOperator2D_Poisson;


// /** @brief a structure for derivatives  in 3D*/
// typedef struct {
//     std::array<std::array<int, 3>, 3> Laplace3D_L  = {{{0,0,0},{0,1,0},{0,0,0}}}; // k-1
//     std::array<std::array<int, 3>, 3> Laplace3D_C  = {{{0,1,0},{1,-6,1},{0,1,0}}}; // k
//     std::array<std::array<int, 3>, 3> Laplace3D_R  = {{{0,0,0},{0,1,0},{0,0,0}}}; // k+1
// } DiscreteOperator3D_Poisson;



/** @brief Discrete Operator Lh */
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

        size_t nb_levels{21}; // maxLevel
        size_t nx_fine{16};
        size_t ny_fine{16};
        size_t nz_fine{16};
        T res_tol{1e-10};
        size_t max_iter{400};
        CycleType cycle_type{multigrid::wCycle};
        size_t npre{1}, npost{1};
        T w{1.0};
        bool fixedNiter{false};
        int Nghost{1};
        BcType bc_condition_type {Periodic};
    };


}


#endif