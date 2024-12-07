#ifndef __MULTIGRID__
#define __MULTIGRID__

#include <vector>
#include <array>
#include "utilities.hpp"
#include "Bcondition.hpp"
#include "ndArray.h"

namespace multigrid
{

    /** @brief restriction operator */
    template <typename T, int Ndim>
    inline void restriction_operator(std::vector<std::vector<T>> &corseData, std::vector<std::vector<T>> &fineData)
    {
    }

    /** @brief interpolation operator */
    template <typename T, int Ndim>
    inline void interpolation_operator(std::vector<std::vector<T>> &corseData, std::vector<std::vector<T>> &fineData)
    {
    }

    /** @brief linear restriction operator */
    template <typename T, int Ndim>
    inline void linear_restriction_operator(ndArray<T, Ndim> &corseData, ndArray<T, Ndim> &fineData)
    {
        size_t nxc,nyc,nzc,nxf,nyf,nzf;
        nxc = corseData.size(0); nyc = corseData.size(1); nzc = (Ndim==3) ? (corseData.size(2)) : 0;
        nxf = fineData.size(0); nyf = fineData.size(1); nzf = (Ndim == 3) ? (fineData.size(2)) : 0;

        if(Ndim == 2)
        {
            // update iner cells
            for(int j = 1; j < nyc-1; j++)
            {
                auto fj = 2*j;
                for(int i = 1; i < nxc - 1; i++)
                {
                    auto fi = 2*i;
                    corseData({i,j}) = (1/64) * (9*fineData({fi, fj}) + 9*fineData({fi + 1, fj}) + 9*fineData({fi + 1, fj + 1}) + 9*fineData({fi , fj + 1})\
                    + 3*fineData({fi - 1, fj}) + 3*fineData({fi + 2, fj }) + 3*fineData({fi + 2, fj + 1})  + 3*fineData({fi - 1 , fj + 1 })\
                    + fineData({fi - 1, fj - 1}) + fineData({fi + 2, fj - 1}) + fineData({fi + 2, fj + 2}) + fineData({fi -1, fj + 2})\
                    + 3*fineData({fi, fj -1 }) + 3*fineData({fi + 1  , 2*j -1}) + 3*fineData({fi + 1, fj + 2}) + 3*fineData({fi , fj + 2})
                    );
                }
            }

            // update boundaries cells (deepending on the BC types e.g Dirichlet, Neumann, Robin, etc..)
        }

        else if(Ndim == 3)
        {
            // update iner cells
            for(int k = 1; k < nzc - 1; k++)
           { 
            auto fk =  2*k;
                for(int j = 1; j < nyc-1; j++)
                {
                    auto fj = 2*j;
                    for(int i = 1; i < nxc - 1; i++)
                    {
                        auto fi = 2*i;
                       corseData({i,j,k}) = (1/512) * (fineData({fi-1,fj+2,fk-1}) + 3*fineData({fi-1,fj+1,fk-1}) + 3*fineData({fi-1,fj,fk-1}) + fineData({fi-1,fj-1,fk-1})\
                                                    + 3*fineData({fi,fj+2,fk-1}) + 9*fineData({fi, fj+1,fk-1}) + 9*fineData({fi,fj,fk-1}) + 3*fineData({fi,fj-1,fk-1}) \
                                                    + 3*fineData({fi+1,fj+2,fk-1}) + 9*fineData({fi+1,fj+1,fk-1}) + 9*fineData({fi+1,fj,fk-1}) + 3*fineData({fi+1,fj-1,fk-1})\
                                                    + fineData({fi+2, fj+2, fk-1}) + 3*fineData({fi+2,fj+1,fk-1}) + 3*fineData({fi+2,fj,fk-1}) + fineData({fi+2,fj-1,fk-1}) \

                                                    + 3*fineData({fi-1,fj+2,fk}) + 9*fineData({fi-1,fj+1,fk}) + 9*fineData({fi-1,fj,fk}) + 3*fineData({fi-1,fj-1,fk}) \
                                                    + 9*fineData({fi,fj+2,fk}) + 27*fineData({fi,fj+1,fk}) + 27*fineData({fi,fj,fk}) + 9*fineData({fi,fj-1,fk}) \
                                                    + 9*fineData({fi+1,fj+2,fk}) + 27*fineData({fi+1,fj+1,fk}) + 27*fineData({fi+1,fj,fk}) + 9*fineData({fi+1,fj-1,fk})\
                                                    + 3*fineData({fi+2,fj+2,fk}) + 9*fineData({fi+2,fj+1,fk}) + 9*fineData({fi+2,fj,fk}) + 3*fineData({fi+2,fj-1,fk})\
                                                    
                                                    + 3*fineData({fi-1,fj+2,fk+1}) + 9*fineData({fi-1,fj+1,fk+1}) + 9*fineData({fi-1,fj,fk+1}) + 3*fineData({fi-1,fj-1,fk+1}) \
                                                    + 9*fineData({fi,fj+2,fk+1}) + 27*fineData({fi,fj+1,fk+1}) + 27*fineData({fi,fj,fk+1}) + 9*fineData({fi,fj-1,fk+1}) \
                                                    + 9*fineData({fi+1,fj+2,fk+1}) + 27*fineData({fi+1,fj+1,fk+1}) + 27*fineData({fi+1,fj,fk+1}) + 9*fineData({fi+1,fj-1,fk+1})\
                                                    + 3*fineData({fi+2,fj+2,fk+1}) + 9*fineData({fi+2,fj+1,fk+1}) + 9*fineData({fi+2,fj,fk+1}) + 3*fineData({fi+2,fj-1,fk+1})\

                                                    + fineData({fi-1,fj+2,fk+2}) + 3*fineData({fi-1,fj+1,fk+2}) + 3*fineData({fi-1,fj,fk+2}) + fineData({fi-1,fj-1,fk+2})\
                                                    + 3*fineData({fi,fj+2,fk+2}) + 9*fineData({fi, fj+1,fk+2}) + 9*fineData({fi,fj,fk+2}) + 3*fineData({fi,fj-1,fk+2}) \
                                                    + 3*fineData({fi+1,fj+2,fk+2}) + 9*fineData({fi+1,fj+1,fk+2}) + 9*fineData({fi+1,fj,fk+2}) + 3*fineData({fi+1,fj-1,fk+2})\
                                                    + fineData({fi+2, fj+2, fk+2}) + 3*fineData({fi+2,fj+1,fk+2}) + 3*fineData({fi+2,fj,fk+2}) + fineData({fi+2,fj-1,fk+2})
                                                    );
                    }
                }
            }

            // updates boundaries cells
        }

    }


        /** @brief khalil restriction operator */
    template <typename T, int Ndim>
    inline void Khalil_restriction_operator(ndArray<T, Ndim> &corseData, ndArray<T, Ndim> &fineData)
    {
        size_t nxc,nyc,nzc,nxf,nyf,nzf;
        nxc = corseData.size(0); nyc = corseData.size(1); nzc = (Ndim==3) ? (corseData.size(2)) : 0;
        nxf = fineData.size(0); nyf = fineData.size(1); nzf = (Ndim == 3) ? (fineData.size(2)) : 0;

        if(Ndim == 2)
        {
            // update iner cells
            for(int j = 1; j < nyc-1; j++)
            {
                auto fj = 2*j;
                for(int i = 1; i < nxc - 1; i++)
                {
                    auto fi = 2*i;
                    corseData({i,j}) = (1/16) * (2*fineData({fi, fj}) + 3*fineData({fi + 1, fj}) + 3*fineData({fi, fj + 1}) + 3*fineData({fi + 1, fj + 1})\
                    + fineData({fi, fj + 2}) + 0*fineData({fi , fj -1}) + 0*fineData({fi + 1, fj + 2}) + fineData({fi + 1, fj - 1}) \
                    + 0*fineData({fi - 1, fj}) + fineData({fi - 1, fj + 1}) + 0*fineData({fi + 2, fj + 1}) + 1*fineData({fi + 2, fj })\
                    + fineData({fi - 1, fj + 2}) + 0*fineData({fi - 1, fj - 1}) + 0*fineData({fi + 2, fj + 2}) + fineData({fi + 1, fj - 1}));
                }
            }

            // update boundaries cells (deepending on the BC types e.g Dirichlet, Neumann, Robin, etc..)
        }

        else if(Ndim == 3)
        {
            // update iner cells
            for(int k = 1; k < nzc - 1; k++)
           { 
            auto fk = 2*k;
                for(int j = 1; j < nyc-1; j++)
                {
                    auto fj = 2*j;
                    for(int i = 1; i < nxc - 1; i++)
                    {
                        auto fi = 2*i;
                       corseData({i,j,k}) = (1/32) * (fineData({fi+1,fj,fk-1}) + fineData({fi+1,fj-1,fk-1}) \
                                                    + fineData({fi+2,fj,fk-1}) + fineData({fi+2,fj-1,fk-1}) \
                                                    + 2*fineData({fi,fj+1,fk}) + 2*fineData({fi,fj,fk}) \
                                                    + 2*fineData({fi+1,fj+1,fk}) + 3*fineData({fi+1,fj,fk}) + fineData({fi+1,fj-1,fk})\
                                                    + fineData({fi+2,fj,fk}) + fineData({fi+2,fj-1,fk})\
                                                    + fineData({fi-1,fj+2,fk+1}) + fineData({fi-1,fj+1,fk+1})  \
                                                    + fineData({fi,fj+2,fk+1}) + 3*fineData({fi,fj+1,fk+1}) + 2*fineData({fi,fj,fk+1}) \
                                                    +  2*fineData({fi+1,fj+1,fk+1}) + 2*fineData({fi+1,fj,fk+1})\
                                                    + fineData({fi-1,fj+2,fk+2}) + fineData({fi-1,fj+1,fk+2})\
                                                    + fineData({fi,fj+2,fk+2}) + fineData({fi, fj+1,fk+2}));
                    }
                }
            }

            // updates boundaries cells
        }

    }


     /** @brief kwak restriction operator */
    template <typename T, int Ndim>
    inline void kwak_restriction_operator(ndArray<T, Ndim> &corseData, ndArray<T, Ndim> &fineData)
    {
        size_t nxc,nyc,nzc,nxf,nyf,nzf;
        nxc = corseData.size(0); nyc = corseData.size(1); nzc = (Ndim==3) ? (corseData.size(2)) : 0;
        nxf = fineData.size(0); nyf = fineData.size(1); nzf = (Ndim == 3) ? (fineData.size(2)) : 0;

        if(Ndim == 2)
        {
            // update iner cells
            for(int j = 1; j < nyc-1; j++)
            {
                auto fj = 2*j;
                for(int i = 1; i < nxc - 1; i++)
                {
                    auto fi = 2*i;
                    corseData({i,j}) = (1/16) * (2*fineData({fi, fj}) + 2*fineData({fi + 1, fj}) + 2*fineData({fi, fj + 1}) + 2*fineData({fi + 1, fj + 1})\
                    + fineData({fi, fj + 2}) + fineData({fi , fj -1}) + fineData({fi + 1, fj + 2}) + fineData({fi + 1, fj - 1}) \
                    + fineData({fi - 1, fj}) + fineData({fi - 1, fj + 1}) + fineData({fi + 2, fj + 1}) + fineData({fi + 2, fj })\
                    + 0*fineData({fi - 1, fj + 2}) + 0*fineData({fi - 1, fj - 1}) + 0*fineData({fi + 2, fj + 2}) + 0*fineData({fi + 1, fj - 1}));
                }
            }

            // update boundaries cells (deepending on the BC types e.g Dirichlet, Neumann, Robin, etc..)
        }

        else if(Ndim == 3)
        {
            // update iner cells
            for(int k = 1; k < nzc - 1; k++)
           { 
            auto fk= 2*k;
                for(int j = 1; j < nyc-1; j++)
                {
                    auto fj = 2*j;
                    for(int i = 1; i < nxc - 1; i++)
                    {
                        auto fi = 2*i;
                        corseData({i,j,k}) = (1/48) * (fineData({i,j+1,k-1}) + fineData({i,j,k-1}) + fineData({i+1,j,k-1}) + fineData({i+1,j+1,k-1}) \

                        + fineData({i-1,j+1,k}) + fineData({i-1,j,k})\
                        + fineData({i,j+2,k}) + 3*fineData({i,j+1,k}) + 3*fineData({i,j,k}) + fineData({i,j-1,k}) \
                        + fineData({i+1,j+2,k}) + 3*fineData({i+1,j+1,k}) + 3*fineData({i+1,j,k}) + fineData({i+1,j-1,k}) \
                        + fineData({i+2,j+1,k}) + fineData({i+2, j,k})\

                        + fineData({i-1,j+1,k+1}) + fineData({i-1,j,k+1})\
                        + fineData({i,j+2,k+1}) + 3*fineData({i,j+1,k+1}) + 3*fineData({i,j,k+1}) + fineData({i,j-1,k+1}) \
                        + fineData({i+1,j+2,k+1}) + 3*fineData({i+1,j+1,k+1}) + 3*fineData({i+1,j,k+1}) + fineData({i+1,j-1,k+1}) \
                        + fineData({i+2,j+1,k+1}) + fineData({i+2, j,k+1})\

                        + fineData({i,j+1,k+2}) + fineData({i,j,k+2}) + fineData({i+1,j,k+2}) + fineData({i+1,j+1,k+2})
                        );
                    }
                }
            }

            // updates boundaries cells
        }

    }



     /** @brief constant restriction operator */
    template <typename T, int Ndim>
    inline void constant_restriction_operator(ndArray<T, Ndim> &corseData, ndArray<T, Ndim> &fineData)
    {
        size_t nxc,nyc,nzc,nxf,nyf,nzf;
        nxc = corseData.size(0); nyc = corseData.size(1); nzc = (Ndim==3) ? (corseData.size(2)) : 0;
        nxf = fineData.size(0); nyf = fineData.size(1); nzf = (Ndim == 3) ? (fineData.size(2)) : 0;

        if(Ndim == 2)
        {
            // update iner cells
            for(int j = 1; j < nyc-1; j++)
            {
                auto fj = 2j;
                for(int i = 1; i < nxc - 1; i++)
                {
                    auto fi = 2*i;
                    corseData({i,j}) = (1/4) * (fineData({fi, fj}) + fineData({fi + 1, fj}) + fineData({fi, fj + 1}) + fineData({fi + 1, fj + 1}));
                }
            }

            // update boundaries cells (deepending on the BC types e.g Dirichlet, Neumann, Robin, etc..)
        }

        else if(Ndim == 3)
        {
            // update iner cells
            for(int k = 1; k < nzc - 1; k++)
           {
            auto fk = 2*k;
                for(int j = 1; j < nyc-1; j++)
                {
                    auto fj = 2*j;
                    for(int i = 1; i < nxc - 1; i++)
                    {
                        auto fi = 2*i;
                        corseData({i,j,k}) = (1/8) * (fineData({fi, fj, fk}) + fineData({fi + 1, fj, fk}) \
                        + fineData({fi , fj + 1, fk}) + fineData({fi, fj, fk + 1}) \
                        + fineData({fi + 1, fj, fk + 1}) + fineData({fi, fj+1, fk+1})\
                        + fineData({fi+1, fj+1, fk}) + fineData({fi + 1, fj+1, fk+1}));
                    }
                }
            }

            // updates boundaries cells
        }

    }

    

    /** @brief constantwise interpolation operator */
    template <typename T, int Ndim>
    inline void constantwise_interpolation_operator(ndArray<T, Ndim> &corseData, ndArray<T, Ndim> &fineData)
    {
        size_t nxc,nyc,nzc,nxf,nyf,nzf;
        nxc = corseData.size(0); nyc = corseData.size(1); nzc = (Ndim==3) ? (corseData.size(2)) : 0;
        nxf = fineData.size(0); nyf = fineData.size(1); nzf = (Ndim == 3) ? (fineData.size(2)) : 0;

        if(Ndim == 2)
        {
            // update iner cells
            for(int j = 0; j < nyc; j++)
            {
                auto  fj = 2*j;
                for(int i = 0; i < nxc; i++)
                {
                    auto fi = 2*i;
                    auto tmp = corseData({i,j});
                    fineData({fi,fj}) = tmp;
                    fineData({fi + 1, fj}) = tmp;
                    fineData({fi + 1,fj + 1}) = tmp;
                    fineData({fi , fj + 1}) = tmp;
                    
                }
            }
        }

        if(Ndim == 3)
        {
            // update iner cells
            for(int k = 0; k < nzc; k++)
           { 
            auto fk = 2*k;
            for(int j = 0; j < nyc; j++)
            {
                auto fj = 2*j;
                for(int i = 0; i < nxc; i++)
                {
                    auto fi = 2*i;
                    auto tmp = corseData({i,j,k});
                    fineData({fi,fj, fk}) = tmp;
                    fineData({fi + 1, fj, fk}) = tmp;
                    fineData({fi + 1,fj + 1, fk}) = tmp;
                    fineData({fi , fj + 1, fk}) = tmp;

                    fineData({fi, fj, fk + 1}) = tmp;
                    fineData({fi + 1, fj, fk + 1}) = tmp;
                    fineData({fi + 1,fj + 1, fk + 1}) = tmp;
                    fineData({fi , fj + 1, fk + 1}) = tmp;
                    
                }
            }
            }
        }


    }



    /** @brief linear interpolation operator */
    // TODO make sure ndArray is initialize to zero
    template <typename T, int Ndim>
    inline void linear_interpolation_operator(ndArray<T, Ndim> &corseData, ndArray<T, Ndim> &fineData)
    {
        size_t nxc,nyc,nzc,nxf,nyf,nzf;
        nxc = corseData.size(0); nyc = corseData.size(1); nzc = (Ndim==3) ? (corseData.size(2)) : 0;
        nxf = fineData.size(0); nyf = fineData.size(1); nzf = (Ndim == 3) ? (fineData.size(2)) : 0;

        if(Ndim == 2)
        {
            // update iner cells
            for(int j = 1; j < nyc-1; j++)
            {
                // for(int i = 1; i < nxc-1; i++)
                // {
                //     auto tmp = (1/16)*corseData({i,j});
                //     fineData({2*i,2*j}) += 9*tmp;
                //     fineData({2*i + 1, 2*j}) += 9*tmp;
                //     fineData({2*i + 1, 2*j +1}) += 9*tmp;
                //     fineData({2*i,2*j + 1}) += 9*tmp;

                //     fineData({2*i-1,2*j}) += 3*tmp;
                //     fineData({2*i+2,2*j}) += 3*tmp;
                //     fineData({2*i+2,2*j+1}) += 3*tmp;
                //     fineData({2*i-1,2*j+1}) += 3*tmp;

                                       
                //     fineData({2*i-1,2*j-1}) += tmp;
                //     fineData({2*i+2,2*j-1}) += tmp;
                //     fineData({2*i+2,2*j+2}) += tmp;
                //     fineData({2*i-1,2*j+2}) += tmp;

                //     fineData({2*i, 2*j-1}) += 3*tmp;
                //     fineData({2*i + 1,2*j -1}) += 3*tmp;
                //     fineData({2*i+1,2*j +2}) += 3*tmp;
                //     fineData({2*i,2*j+2}) += 3*tmp;
                // }
                auto fj = 2*j;
                for(int i = 1; i < nxc-1; i++)
                {
                    auto fi = 2*i;
                    fineData({fi,fj})  = (1/16) * (9*corseData({i,j})  + 3 * corseData({i-1,j}) + 3 * corseData({i,j-1}) + corseData({i-1,j-1}));
                    fineData({fi+1,fj}) = (1/16) * (9*corseData({i,j}) + 3 * corseData({i+1,j}) + 3 * corseData({i,j-1}) + corseData({i+1,j-1}));
                    fineData({fi,fj+1}) = (1/16) * (9*corseData({i,j}) + 3 * corseData({i,j+1}) + 3 * corseData({i-1,j}) + corseData({i-1,j+1}));
                    fineData({fi+1,fj+1}) = (1/16) * (9*corseData({i,j}) + 3 * corseData({i,j+1}) + 3 * corseData({i+1, j}) + corseData({i+1,j+1}));
                }
            }

            // TODO BC to add on
        }

        if(Ndim == 3)
        {
            for(int k = 1; k < nzc-1; k++)
            {

                auto fk = 2*k;
                for(int j = 1; j < nyc-1; j++)
                {
                    auto fj = 2*j;
                    for(int i = 1; i < nxc-1; i++)
                    {
                        auto fi = 2*i;





                        fineData({fi,fj,fk}) = (1/64) * (27*corseData({i,j,k}) + 9 *corseData({i-1,j,k}) + 3*corseData({i-1,j-1,k}) + 9*corseData({i,j-1,k}) \
                        + 3 * corseData({i-1,j,k-1}) + corseData({i-1,j-1,k-1}) + 3 * corseData({i,j-1,k-1}) + 9 * corseData({i,j,k-1}));

                        fineData({fi+1,fj,fk}) = (1/64) *( 27*corseData(i,j,k)  + 9*corseData({i,j-1,k}) + 3*corseData({i+1,j-1,k}) + 9 * corseData(i+1,j,k) \
                        + 3 * corseData({i,j-1,k-1}) + corseData({i+1,j-1,k-1}) + 9 * corseData({i,j,k-1}) + 3 * corseData({i+1,j,k-1}));

                        fineData({fi,fj+1,fk}) = (1/64) * (27 * corseData(i,j,k) + 9 *corseData({i,j+1,k}) + 9*corseData({i-1,j,k}) + 3 * corseData({i-1,j+1,k})\
                        + 9 * corseData({i,j,k-1}) + 3 * corseData({i,j+1,k-1}) + 3 * corseData({i-1,j,k-1}) + corseData({i-1,j+1,k-1}));

                        fineData({fi+1,fj+1, fk}) = (1/64) * ( 27 * corseData(i,j,k) + 9 * corseData({i+1,j,k}) + 9 * corseData({i,j+1,k}) + 3 * corseData(i+1,j+1,k) \
                        + 9 * corseData({i,j,k-1}) + 3 * corseData({i+1,j,k-1}) + 3 * corseData({i,j+1,k-1}) + corseData(i+1, j+1, k-1));

                        fineData({fi,fj,fk+1}) = (1/64) * ( 27 * corseData({i,j,k}) + 9*corseData({i,j-1,k}) + 9 * corseData({i-1,j,k}) + 3 * corseData({i-1,j-1,k}) \
                        + 9 * corseData({i,j,k+1}) + 3 * corseData({i,j-1,k+1}) + 3*corseData({i-1,j,k+1}) + corseData({i-1,j-1,k+1}));

                        fineData({fi+1,fj,fk+1}) = (1/64) * ( 27 * corseData({i,j,k}) + 9 * corseData({i+1,j,k}) + 9 * corseData({i,j-1,k}) + 3 * corseData(i+1,j-1,k)\
                        + 9 * corseData({i,j,k+1}) + 3 * corseData({i+1,j,k+1}) + 3 * corseData(i,j-1,k+1) + corseData({i+1,j-1,k+1}));

                        fineData({fi,fj+1,fk+1}) = (1/64) * ( 27 * corseData({i,j,k}) + 9 * corseData({i,j+1,k}) + 9 * corseData(i-1,j,k) + 3 * corseData(i-1,j+1,k) \
                        + 9 * corseData({i,j,k+1}) + 3 * corseData({i-1,j,k+1}) + 3 * corseData({i,j+1,k+1}) + corseData({i-1, j+1, k+1}));

                        fineData({fi+1,fj+1,fk+1}) = (1/64) * (27 * corseData({i,j,k}) + 9 * corseData({i+1,j,k}) + 9*corseData({i,j+1,k}) + corseData({i+1,j+1,k}) \
                        + 9 * corseData({i,j,k+1}) + 3 * corseData({i,j+1,k+1}) + 3*corseData({i+1,j,k+1}) + corseData({i+1,j+1,k+1}));

                    // 

                    //     auto tmp = (1/64)*corseData({i,j,k });
                    // {  
                    //     fineData({2*i,2*j, 2*k - 1}) += 9*tmp;
                    //     fineData({2*i + 1, 2*j,  2*k - 1}) += 9*tmp;
                    //     fineData({2*i + 1, 2*j +1,  2*k - 1}) += 9*tmp;
                    //     fineData({2*i,2*j + 1,  2*k - 1}) += 9*tmp;

                    //     fineData({2*i-1,2*j,  2*k - 1}) += 3*tmp;
                    //     fineData({2*i+2,2*j,  2*k - 1}) += 3*tmp;
                    //     fineData({2*i+2,2*j+1,  2*k - 1}) += 3*tmp;
                    //     fineData({2*i-1,2*j+1,  2*k - 1}) += 3*tmp;

                                        
                    //     fineData({2*i-1,2*j-1,  2*k - 1}) += tmp;
                    //     fineData({2*i+2,2*j-1,  2*k - 1}) += tmp;
                    //     fineData({2*i+2,2*j+2,  2*k - 1}) += tmp;
                    //     fineData({2*i-1,2*j+2 , 2*k - 1}) += tmp;

                    //     fineData({2*i, 2*j-1,  2*k - 1}) += 3*tmp;
                    //     fineData({2*i + 1,2*j -1,  2*k - 1}) += 3*tmp;
                    //     fineData({2*i+1,2*j +2,  2*k - 1}) += 3*tmp;
                    //     fineData({2*i,2*j+2,  2*k - 1}) += 3*tmp;
                        
                    // }

            
                    // {  
                    //     fineData({2*i,2*j,  2*k}) +=        27*tmp;
                    //     fineData({2*i + 1, 2*j,  2*k }) +=  27*tmp;
                    //     fineData({2*i + 1, 2*j +1,  2*k})+= 27*tmp;
                    //     fineData({2*i,2*j + 1,  2*k}) +=    27*tmp;

                    //     fineData({2*i-1,2*j,  2*k})  += 9*tmp;
                    //     fineData({2*i+2,2*j,  2*k})  += 9*tmp;
                    //     fineData(2*i+2,2*j+1,  2*k ) += 9*tmp;
                    //     fineData(2*i-1,2*j+1,  2*k ) += 9*tmp;

                                        
                    //     fineData({2*i-1,2*j-1,  2*k }) += 3*tmp;
                    //     fineData({2*i+2,2*j-1,  2*k }) += 3*tmp;
                    //     fineData({2*i+2,2*j+2,  2*k }) += 3*tmp;
                    //     fineData({2*i-1,2*j+2,  2*k }) += 3*tmp;

                    //     fineData({2*i, 2*j-1,      2*k}) += 9*tmp;
                    //     fineData({2*i + 1,2*j -1,  2*k}) += 9*tmp;
                    //     fineData({2*i+1,2*j +2,    2*k}) += 9*tmp;
                    //     fineData({2*i,2*j+2,       2*k}) += 9*tmp;
                        
                    // }

                       
                    // {  
                    //     fineData({2*i,     2*j,     2*k + 1}) += 27*tmp;
                    //     fineData({2*i + 1, 2*j,     2*k + 1}) += 27*tmp;
                    //     fineData({2*i + 1, 2*j +1,  2*k + 1}) += 27*tmp;
                    //     fineData({2*i,     2*j + 1, 2*k + 1}) += 27*tmp;

                    //     fineData({2*i-1,2*j, 2*k + 1})   += 9*tmp;
                    //     fineData({2*i+2,2*j, 2*k + 1})   += 9*tmp;
                    //     fineData({2*i+2,2*j+1, 2*k + 1}) += 9*tmp;
                    //     fineData({2*i-1,2*j+1, 2*k + 1}) += 9*tmp;

                                        
                    //     fineData({2*i-1,2*j-1, 2*k + 1}) += 3*tmp;
                    //     fineData({2*i+2,2*j-1, 2*k + 1}) += 3*tmp;
                    //     fineData({2*i+2,2*j+2, 2*k + 1}) += 3*tmp;
                    //     fineData({2*i-1,2*j+2, 2*k + 1}) += 3*tmp;

                    //     fineData({2*i, 2*j-1,     2*k + 1}) += 9*tmp;
                    //     fineData({2*i + 1,2*j -1, 2*k + 1}) += 9*tmp;
                    //     fineData({2*i+1,2*j +2,   2*k + 1}) += 9*tmp;
                    //     fineData({2*i,2*j+2,      2*k + 1}) += 9*tmp;
                        
                    // }
                    

                    // {  
                    //     fineData({2*i,2*j, 2*k + 2})         += 9*tmp;
                    //     fineData({2*i + 1, 2*j, 2*k + 2})    += 9*tmp;
                    //     fineData({2*i + 1, 2*j +1, 2*k + 2}) += 9*tmp;
                    //     fineData({2*i,2*j + 1, 2*k + 2})     += 9*tmp;

                    //     fineData({2*i-1,2*j, 2*k + 2})   += 3*tmp;
                    //     fineData({2*i+2,2*j, 2*k + 2})   += 3*tmp;
                    //     fineData({2*i+2,2*j+1, 2*k + 2}) += 3*tmp;
                    //     fineData({2*i-1,2*j+1, 2*k + 2}) += 3*tmp;

                                        
                    //     fineData({2*i-1,2*j-1, 2*k + 2}) += tmp;
                    //     fineData({2*i+2,2*j-1, 2*k + 2}) += tmp;
                    //     fineData({2*i+2,2*j+2, 2*k + 2}) += tmp;
                    //     fineData({2*i-1,2*j+2, 2*k + 2}) += tmp;

                    //     fineData({2*i, 2*j-1, 2*k + 2})     += 3*tmp;
                    //     fineData({2*i + 1,2*j -1, 2*k + 2}) += 3*tmp;
                    //     fineData({2*i+1,2*j +2, 2*k + 2})   += 3*tmp;
                    //     fineData({2*i,2*j+2, 2*k + 2})      += 3*tmp;
                        
                    // }

                    // 
                    
                    }
                }
            }

            // TODO : add BC
            
        }
    }

    /** @brief Khalil interpolation operator */
    // TODO make sure ndArray is initialize to zero
    template <typename T, int Ndim>
    inline void khalil_interpolation_operator(ndArray<T, Ndim> &corseData, ndArray<T, Ndim> &fineData)
    {
        size_t nxc,nyc,nzc,nxf,nyf,nzf;
        nxc = corseData.size(0); nyc = corseData.size(1); nzc = (Ndim==3) ? (corseData.size(2)) : 0;
        nxf = fineData.size(0); nyf = fineData.size(1); nzf = (Ndim == 3) ? (fineData.size(2)) : 0;

        if(Ndim == 2)
        {
            // update iner cells
            for(int j = 1; j < nyc-1; j++)
            {
                // for(int i = 1; i < nxc-1; i++)
                // {
                //     auto tmp = (1/4)*corseData({i,j});
                //     fineData({2*i,2*j}) += 2*tmp;
                //     fineData({2*i + 1, 2*j}) += 3*tmp;
                //     fineData({2*i + 1, 2*j +1}) += 2*tmp;
                //     fineData({2*i,2*j + 1}) += 3*tmp;

                //     fineData({2*i-1,2*j}) += 0*tmp;
                //     fineData({2*i+2,2*j}) += tmp;
                //     fineData({2*i+2,2*j+1}) += 0*tmp;
                //     fineData({2*i-1,2*j+1}) += tmp;

                                       
                //     fineData({2*i-1,2*j-1}) += 0*tmp;
                //     fineData({2*i+2,2*j-1}) += tmp;
                //     fineData({2*i+2,2*j+2}) += 0*tmp;
                //     fineData({2*i-1,2*j+2}) += tmp;

                //     fineData({2*i, 2*j-1}) += 0*tmp;
                //     fineData({2*i + 1,2*j -1}) += tmp;
                //     fineData({2*i+1,2*j +2}) +=0*tmp;
                //     fineData({2*i,2*j+2}) += tmp;
                   
                // }
                auto fj = j*2;

                for(int i = 1; i < nxc-1; i++)
                {
                    auto fi = 2*i;
                    fineData({fi,fj}) = (1/4) * (2*corseData({i,j}));
                    fineData({fi+1,fj}) = (1/4) * (3*corseData({i,j}) + corseData({i+1,j}) + corseData({i,j-1}) + corseData({i+1,j-1}));
                    fineData({fi,fj+1}) = (1/4) * (3*corseData({i,j}) + corseData({i-1,j}) + corseData({i,j+1}) + corseData({i-1,j+1}));
                    fineData({fi+1,fj+1}) = (1/4) * (2*corseData(i,j)) ;
                   
                }
            }

            // TODO BC to add on
        }

        if(Ndim == 3)
        {
            for(int k = 1; k < nzc-1; k++)
            {
                auto fk = 2 * k;
                for(int j = 1; j < nyc-1; j++)
                {
                    auto fj = 2 * j;
                    for(int i = 1; i < nxc-1; i++)
                    {
                        auto fi = 2 * i;


                        fineData({fi,fj,fk}) = (1/4) * (2*corseData({i,j,k}));

                        fineData({fi+1,fj,fk}) = (1/4) *( 3*corseData(i,j,k)  + corseData({i,j-1,k}) + corseData({i+1,j-1,k}) + corseData(i+1,j,k) \
                        + corseData({i,j-1,k-1}) + corseData({i+1,j-1,k-1}) + corseData({i,j,k-1}) + corseData({i+1,j,k-1}));

                        fineData({fi,fj+1,fk}) = (1/4) * (2 * corseData(i,j,k));

                        fineData({fi+1,fj+1, fk}) = (1/4) * (2 * corseData(i,j,k));

                        fineData({fi,fj,fk+1}) = (1/4) * (2 * corseData({i,j,k}));

                        fineData({fi+1,fj,fk+1}) = (1/4) * (2 * corseData({i,j,k}));

                        fineData({fi,fj+1,fk+1}) = (1/4) * ( 3 * corseData({i,j,k}) +  corseData({i,j+1,k}) + corseData(i-1,j,k) + corseData(i-1,j+1,k) \
                        +  corseData({i,j,k+1}) + corseData({i-1,j,k+1}) + corseData({i,j+1,k+1}) + corseData({i-1, j+1, k+1}));

                        fineData({fi+1,fj+1,fk+1}) = (1/4) * (2 * corseData({i,j,k}));



                    //     auto tmp = (1/4)*corseData({i,j,k });
                    // {  
                    //     fineData({2*i,2*j, 2*k - 1}) +=          0*tmp;
                    //     fineData({2*i + 1, 2*j,  2*k - 1})  +=   tmp;
                    //     fineData({2*i + 1, 2*j +1,  2*k - 1}) += 0*tmp;
                    //     fineData({2*i,2*j + 1,  2*k - 1}) +=     0*tmp;

                    //     fineData({2*i-1,2*j,  2*k - 1}) +=   0*tmp;
                    //     fineData({2*i+2,2*j,  2*k - 1}) +=   tmp;
                    //     fineData({2*i+2,2*j+1,  2*k - 1}) += 0*tmp;
                    //     fineData({2*i-1,2*j+1,  2*k - 1}) += 0*tmp;

                                        
                    //     fineData({2*i-1,2*j-1,  2*k - 1}) += 0*tmp;
                    //     fineData({2*i+2,2*j-1,  2*k - 1}) += tmp;
                    //     fineData({2*i+2,2*j+2,  2*k - 1}) += 0*tmp;
                    //     fineData({2*i-1,2*j+2 , 2*k - 1}) += 0*tmp;

                    //     fineData({2*i, 2*j-1,  2*k - 1}) +=     0*tmp;
                    //     fineData({2*i + 1,2*j -1,  2*k - 1}) += tmp;
                    //     fineData({2*i+1,2*j +2,  2*k - 1}) +=   0*tmp;
                    //     fineData({2*i,2*j+2,  2*k - 1}) +=      0*tmp;
                        
                    // }

            
                    // {  
                    //     fineData({2*i,2*j,  2*k}) +=        2*tmp;
                    //     fineData({2*i + 1, 2*j,  2*k }) +=  3*tmp;
                    //     fineData({2*i + 1, 2*j +1,  2*k})+= 2*tmp;
                    //     fineData({2*i,2*j + 1,  2*k}) +=    2*tmp;

                    //     fineData({2*i-1,2*j,  2*k})  += 0*tmp;
                    //     fineData({2*i+2,2*j,  2*k})  += tmp;
                    //     fineData(2*i+2,2*j+1,  2*k ) += 0*tmp;
                    //     fineData(2*i-1,2*j+1,  2*k ) += 0*tmp;

                                        
                    //     fineData({2*i-1,2*j-1,  2*k }) += 0*tmp;
                    //     fineData({2*i+2,2*j-1,  2*k }) += tmp;
                    //     fineData({2*i+2,2*j+2,  2*k }) += 0*tmp;
                    //     fineData({2*i-1,2*j+2,  2*k }) += 0*tmp;

                    //     fineData({2*i, 2*j-1,      2*k}) += 0*tmp;
                    //     fineData({2*i + 1,2*j -1,  2*k}) += tmp;
                    //     fineData({2*i+1,2*j +2,    2*k}) += 0*tmp;
                    //     fineData({2*i,2*j+2,       2*k}) += 0*tmp;
                        
                    // }

                       
                    // {  
                    //     fineData({2*i,     2*j,     2*k + 1}) += 2*tmp;
                    //     fineData({2*i + 1, 2*j,     2*k + 1}) += 2*tmp;
                    //     fineData({2*i + 1, 2*j +1,  2*k + 1}) += 2*tmp;
                    //     fineData({2*i,     2*j + 1, 2*k + 1}) += 3*tmp;

                    //     fineData({2*i-1,2*j, 2*k + 1})   += 0*tmp;
                    //     fineData({2*i+2,2*j, 2*k + 1})   += 0*tmp;
                    //     fineData({2*i+2,2*j+1, 2*k + 1}) += 0*tmp;
                    //     fineData({2*i-1,2*j+1, 2*k + 1}) += tmp;

                                        
                    //     fineData({2*i-1,2*j-1, 2*k + 1}) += 0*tmp;
                    //     fineData({2*i+2,2*j-1, 2*k + 1}) += 0*tmp;
                    //     fineData({2*i+2,2*j+2, 2*k + 1}) += 0*tmp;
                    //     fineData({2*i-1,2*j+2, 2*k + 1}) += tmp;

                    //     fineData({2*i, 2*j-1,     2*k + 1}) += 0*tmp;
                    //     fineData({2*i + 1,2*j -1, 2*k + 1}) += 0*tmp;
                    //     fineData({2*i+1,2*j +2,   2*k + 1}) += 0*tmp;
                    //     fineData({2*i,2*j+2,      2*k + 1}) += tmp;
                        
                    // }
                    

                    // {  
                    //     fineData({2*i,2*j, 2*k + 2})         += 0*tmp;
                    //     fineData({2*i + 1, 2*j, 2*k + 2})    += 0*tmp;
                    //     fineData({2*i + 1, 2*j +1, 2*k + 2}) += 0*tmp;
                    //     fineData({2*i,2*j + 1, 2*k + 2})     += tmp;

                    //     fineData({2*i-1,2*j, 2*k + 2})   += 0*tmp;
                    //     fineData({2*i+2,2*j, 2*k + 2})   += 0*tmp;
                    //     fineData({2*i+2,2*j+1, 2*k + 2}) += 0*tmp;
                    //     fineData({2*i-1,2*j+1, 2*k + 2}) += tmp;

                                        
                    //     fineData({2*i-1,2*j-1, 2*k + 2}) += 0*tmp;
                    //     fineData({2*i+2,2*j-1, 2*k + 2}) += 0*tmp;
                    //     fineData({2*i+2,2*j+2, 2*k + 2}) += 0*tmp;
                    //     fineData({2*i-1,2*j+2, 2*k + 2}) += tmp;

                    //     fineData({2*i, 2*j-1, 2*k + 2})     += 0*tmp;
                    //     fineData({2*i + 1,2*j -1, 2*k + 2}) += 0*tmp;
                    //     fineData({2*i+1,2*j +2, 2*k + 2})   += 0*tmp;
                    //     fineData({2*i,2*j+2, 2*k + 2})      += tmp;
                        
                    // }
                    
                    }
                }
            }

            // TODO : add BC
        }
    }


    /** @brief Kwak interpolation operator */
    // TODO make sure ndArray is initialize to zero
    template <typename T, int Ndim>
    inline void kwak_interpolation_operator(ndArray<T, Ndim> &corseData, ndArray<T, Ndim> &fineData)
    {
        size_t nxc,nyc,nzc,nxf,nyf,nzf;
        nxc = corseData.size(0); nyc = corseData.size(1); nzc = (Ndim==3) ? (corseData.size(2)) : 0;
        nxf = fineData.size(0); nyf = fineData.size(1); nzf = (Ndim == 3) ? (fineData.size(2)) : 0;

        if(Ndim == 2)
        {
            // update iner cells
            for(int j = 1; j < nyc-1; j++)
            {
        
                // for(int i = 1; i < nxc-1; i++)
                // {
                //     auto tmp = (1/4)*corseData({i,j});
                //     fineData({2*i,2*j}) += 2*tmp;
                //     fineData({2*i + 1, 2*j}) += 2*tmp;
                //     fineData({2*i + 1, 2*j +1}) += 2*tmp;
                //     fineData({2*i,2*j + 1}) += 2*tmp;

                //     fineData({2*i-1,2*j}) += tmp;
                //     fineData({2*i+2,2*j}) += tmp;
                //     fineData({2*i+2,2*j+1}) += tmp;
                //     fineData({2*i-1,2*j+1}) += tmp;

                                       
                //     fineData({2*i-1,2*j-1}) += 0*tmp;
                //     fineData({2*i+2,2*j-1}) += 0*tmp;
                //     fineData({2*i+2,2*j+2}) += 0*tmp;
                //     fineData({2*i-1,2*j+2}) += 0*tmp;

                //     fineData({2*i, 2*j-1}) += tmp;
                //     fineData({2*i + 1,2*j -1}) += tmp;
                //     fineData({2*i+1,2*j +2}) += tmp;
                //     fineData({2*i,2*j+2}) += tmp;
                   
                // }
                auto fj = 2*j;
                for(int i = 1; i < nxc-1; i++)
                {
                    auto fi = 2*i;
                    fineData({fi,fj}) = (1/4) * (2*corseData({i,j}) + corseData({i-1,j}) + corseData({i,j-1}));
                    fineData({fi+1,fj}) = (1/4) * (2*corseData({i,j}) + corseData({i,j-1}) + corseData({i+1,j}));
                    fineData({fi,fj+1}) = (1/4) * (2*corseData({i,j}) + corseData({i-1, j}) + corseData({i,j+1}));
                    fineData({fi+1,fj+1}) = (1/4) * (2*corseData({i,j}) + corseData({i+1,j}) + corseData({i+1,j+1}));
                   
                }
            }

            // TODO BC to add on
        }

        if(Ndim == 3)
        {
            for(int k = 1; k < nzc-1; k++)
            {
                auto fk = 2*k;
                for(int j = 1; j < nyc-1; j++)
                {
                    auto fj = 2*j;
                    for(int i = 1; i < nxc-1; i++)
                    {
                        auto fi = 2*i;

                        fineData({fi,fj,fk}) = (1/6) *(3*corseData({i,j,k}) + corseData({i-1,j,k}) + corseData({i-1,j-1,k}) + corseData({i,j,k-1}));

                        fineData({fi+1,fj,fk}) = (1/6) *(3*corseData({i,j,k})  + corseData({i,j-1,k}) +  corseData({i+1,j,k}) + corseData({i,j,k-1}));

                        fineData({fi,fj+1,fk}) = (1/6) *(3 * corseData({i,j,k}) + corseData({i,j+1,k}) + corseData({i-1,j,k}) + corseData({i,j,k-1}));

                        fineData({fi+1,fj+1, fk}) = (1/6) * (3* corseData(i,j,k) +corseData({i+1,j,k}) +corseData({i,j+1,k}) + corseData({i,j,k-1}));

                        fineData({fi,fj,fk+1}) = (1/6) * (3 * corseData({i,j,k}) + corseData({i,j-1,k}) +  corseData({i-1,j,k}) + corseData({i,j,k+1}));

                        fineData({fi+1,fj,fk+1}) = (1/6) * (3 * corseData({i,j,k}) +  corseData({i+1,j,k}) + corseData({i,j-1,k}) +  corseData({i,j,k+1}));

                        fineData({fi,fj+1,fk+1}) = (1/6) * (3 * corseData({i,j,k}) + corseData({i,j+1,k}) +  corseData(i-1,j,k) + corseData({i,j,k+1}));

                        fineData({fi+1,fj+1,fk+1}) = (1/6) * (3 * corseData({i,j,k}) + corseData({i+1,j,k}) + corseData({i,j+1,k}) + corseData({i,j,k+1}));

                    //     auto tmp = (1/6)*corseData({i,j,k });
                    // {  
                    //     fineData({2*i,2*j, 2*k - 1}) +=          tmp;
                    //     fineData({2*i + 1, 2*j,  2*k - 1}) +=    tmp;
                    //     fineData({2*i + 1, 2*j +1,  2*k - 1}) += tmp;
                    //     fineData({2*i,2*j + 1,  2*k - 1}) +=     tmp;

                    //     fineData({2*i-1,2*j,  2*k - 1}) +=    0*tmp;
                    //     fineData({2*i+2,2*j,  2*k - 1}) +=    0*tmp;
                    //     fineData({2*i+2,2*j+1,  2*k - 1}) +=  0*tmp;
                    //     fineData({2*i-1,2*j+1,  2*k - 1}) +=  0*tmp;

                                        
                    //     fineData({2*i-1,2*j-1,  2*k - 1}) += 0*tmp;
                    //     fineData({2*i+2,2*j-1,  2*k - 1}) += 0*tmp;
                    //     fineData({2*i+2,2*j+2,  2*k - 1}) += 0*tmp;
                    //     fineData({2*i-1,2*j+2 , 2*k - 1}) += 0*tmp;

                    //     fineData({2*i, 2*j-1,  2*k - 1}) +=     0*tmp;
                    //     fineData({2*i + 1,2*j -1,  2*k - 1}) += 0*tmp;
                    //     fineData({2*i+1,2*j +2,  2*k - 1}) +=   0*tmp;
                    //     fineData({2*i,2*j+2,  2*k - 1}) +=      0*tmp;
                        
                    // }

            
                    // {  
                    //     fineData({2*i,2*j,  2*k}) +=        3*tmp;
                    //     fineData({2*i + 1, 2*j,  2*k }) +=  3*tmp;
                    //     fineData({2*i + 1, 2*j +1,  2*k})+= 3*tmp;
                    //     fineData({2*i,2*j + 1,  2*k}) +=    3*tmp;

                    //     fineData({2*i-1,2*j,  2*k})  += tmp;
                    //     fineData({2*i+2,2*j,  2*k})  += tmp;
                    //     fineData(2*i+2,2*j+1,  2*k ) += tmp;
                    //     fineData(2*i-1,2*j+1,  2*k ) += tmp;

                                        
                    //     fineData({2*i-1,2*j-1,  2*k }) += 0*tmp;
                    //     fineData({2*i+2,2*j-1,  2*k }) += 0*tmp;
                    //     fineData({2*i+2,2*j+2,  2*k }) += 0*tmp;
                    //     fineData({2*i-1,2*j+2,  2*k }) += 0*tmp;

                    //     fineData({2*i, 2*j-1,      2*k}) += tmp;
                    //     fineData({2*i + 1,2*j -1,  2*k}) += tmp;
                    //     fineData({2*i+1,2*j +2,    2*k}) += tmp;
                    //     fineData({2*i,2*j+2,       2*k}) += tmp;
                        
                    // }

                       
                    // {  
                    //     fineData({2*i,     2*j,     2*k + 1}) += 3*tmp;
                    //     fineData({2*i + 1, 2*j,     2*k + 1}) += 3*tmp;
                    //     fineData({2*i + 1, 2*j +1,  2*k + 1}) += 3*tmp;
                    //     fineData({2*i,     2*j + 1, 2*k + 1}) += 3*tmp;

                    //     fineData({2*i-1,2*j, 2*k + 1})   += tmp;
                    //     fineData({2*i+2,2*j, 2*k + 1})   += tmp;
                    //     fineData({2*i+2,2*j+1, 2*k + 1}) += tmp;
                    //     fineData({2*i-1,2*j+1, 2*k + 1}) += tmp;

                                        
                    //     fineData({2*i-1,2*j-1, 2*k + 1}) += 0*tmp;
                    //     fineData({2*i+2,2*j-1, 2*k + 1}) += 0*tmp;
                    //     fineData({2*i+2,2*j+2, 2*k + 1}) += 0*tmp;
                    //     fineData({2*i-1,2*j+2, 2*k + 1}) += 0*tmp;

                    //     fineData({2*i, 2*j-1,     2*k + 1}) += tmp;
                    //     fineData({2*i + 1,2*j -1, 2*k + 1}) += tmp;
                    //     fineData({2*i+1,2*j +2,   2*k + 1}) += tmp;
                    //     fineData({2*i,2*j+2,      2*k + 1}) += tmp;
                        
                    // }
                    

                    // {  
                    //     fineData({2*i,2*j, 2*k + 2})         += tmp;
                    //     fineData({2*i + 1, 2*j, 2*k + 2})    += tmp;
                    //     fineData({2*i + 1, 2*j +1, 2*k + 2}) += tmp;
                    //     fineData({2*i,2*j + 1, 2*k + 2})     += tmp;

                    //     fineData({2*i-1,2*j, 2*k + 2})   += 0*tmp;
                    //     fineData({2*i+2,2*j, 2*k + 2})   += 0*tmp;
                    //     fineData({2*i+2,2*j+1, 2*k + 2}) += 0*tmp;
                    //     fineData({2*i-1,2*j+1, 2*k + 2}) += 0*tmp;

                                        
                    //     fineData({2*i-1,2*j-1, 2*k + 2}) += 0*tmp;
                    //     fineData({2*i+2,2*j-1, 2*k + 2}) += 0*tmp;
                    //     fineData({2*i+2,2*j+2, 2*k + 2}) += 0*tmp;
                    //     fineData({2*i-1,2*j+2, 2*k + 2}) += 0*tmp;

                    //     fineData({2*i, 2*j-1, 2*k + 2})     += 0*tmp;
                    //     fineData({2*i + 1,2*j -1, 2*k + 2}) += 0*tmp;
                    //     fineData({2*i+1,2*j +2, 2*k + 2})   += 0*tmp;
                    //     fineData({2*i,2*j+2, 2*k + 2})      += 0*tmp;
                        
                    // }
                    
                    }
                }
            }

            // update BC
        }
    }



    /** @brief structure to set initial parameter for the multigrid class */
    template <typename T>
    struct Settings
    {
        // T aspectRatio{1};
        size_t numberOfGrids{8}; // maxLevel
        size_t minimumResolution{4};
        T residualTolerance{1e-10};
        size_t maxIterations{400};
        CycleType cycleType{multigrid::wCycle};
        size_t Npre{1}, Npost{1}, maxLevel{21}, Nvcyle_cvg{20};
    };

    /** @brief a stack to track grids at all level */
    template <typename T, int Ndim>
    class Stack : public std::vector<ndArray<T, Ndim>>
    {
    public:
        Stack(const Settings<T> &s)
        {
            nGrids = s.numberOfGrids;
            mini_resolution = s.minimumResolution;
            finestLevel = (nGrids - 1);
            this->resize(nGrids);
            // set coarsest grid size
            int nx{mini_resolution}, ny{mini_resolution}, nz{mini_resolution};
            if (Ndim == 2)
                nz = 0;

            // Data set for all levels
            for (auto level = 1; level <= finestLevel; level++)
            {
                nx = 2 * (nx - 1) + 1;
                ny = 2 * (ny - 1) + 1;
                nz = (Ndim == 3) ? (2 * (nz - 1) + 1) : 0;
                auto lin_size = (Ndim == 1) ? (nx * ny) : (nx * ny * nz);
                // (*this)[level] = std::vector<std::vector<T>>(std::vector<T>(0, nx), ny);
                std::vector<T> tmp_init_data(0, lin_size);
                (*this)[level] = ndArray<T, Ndim>(tmp_init_data.data(), {nx, ny, 1}, true);
            }

            // TODO Boundary conditions for all levels
        };
        virtual ~Stack() {};

        size_t finestLevel;                    // finest level
        static const size_t coarsestLevel = 0; // coarsest level

        /** @brief function to coarse a grid from level l to level l-1 */
        inline void coarsen(size_t level);

        /** @brief function to coarse a grid from level l to level l-1 and store result in coarseLevelData */
        inline void coarsen(size_t level, std::vector<std::vector<T>> &coarseLevelData);
        /** @brief function to coarse a grid from level l to level l-1 and store result in coarseLevelData */
        inline void coarsen(size_t level, ndArray<T, Ndim> &coarseLevelData);

        /** @brief function to refine a grid from level l-1 to level l */
        inline void refine(size_t level);
        /** @brief function to coarse a grid from level l to level l-1 and store result in fineLevelData */
        inline void refine(size_t level, std::vector<std::vector<T>> &fineLevelData);
        /** @brief function to coarse a grid from level l to level l-1 and store result in fineLevelData */
        inline void refine(size_t level, ndArray<T, Ndim> &fineLevelData);
        BoundaryConditions bc; // boundary condition

    private:
        const size_t nGrids; // number of grid levels required
        int mini_resolution; // minimum resolution
    };

    template <typename T, int Ndim>
    inline void Stack<T, Ndim>::coarsen(size_t level)
    {
        restriction_operator<T, Ndim>((*this)[level - 1], (*this)[level]);
    }

    template <typename T, int Ndim>
    inline void Stack<T, Ndim>::coarsen(size_t level, std::vector<std::vector<T>> &coarseLevelData)
    {
        restriction_operator(coarseLevelData, (*this)[level]);
    }

    template <typename T, int Ndim>
    inline void Stack<T, Ndim>::coarsen(size_t level, ndArray<T, Ndim> &coarseLevelData)
    {
        restriction_operator(coarseLevelData, (*this)[level]);
    }

    template <typename T, int Ndim>
    inline void Stack<T, Ndim>::refine(size_t level)
    {
        interpolation_operator((*this)[level], (*this)[level + 1]);
    }

    template <typename T, int Ndim>
    inline void Stack<T, Ndim>::refine(size_t level, std::vector<std::vector<T>> &fineLevelData)
    {
        interpolation_operator((*this)[level], fineLevelData);
    }

    template <typename T, int Ndim>
    inline void Stack<T, Ndim>::refine(size_t level, ndArray<T, Ndim> &fineLevelData)
    {
        interpolation_operator((*this)[level], fineLevelData);
    }

    template <typename T, int Ndim>
    class MultigridBase
    {
    public:
        MultigridBase(const Settings<T> &settings);
        virtual ~MultigridBase() {};

        /** @brief get the output of the multigrid method */
        inline ndArray<T, Ndim> &get_result();
        // inline std::vector<std::vector<T>> &get_result();

        /** @brief set the initial guess (e.g potential in the case of Poisson equation) */
        template <typename U>
        inline void initial_guess(U arg);
        // inline std::vector<std::vector<T>> &source_term();

        /** @brief get the soruce term */
        inline ndArray<T, Ndim> &source_term();

        /** @brief set the source term with the given argument */
        template <typename U>
        inline void source_term(U arg);

        /** @brief Compute L_{H}U_{H} on a given grid of resolution H */
        inline void evaluate_operator(size_t level, std::vector<std::vector<T>> &result);

        /** @brief Compute compute the residual */
        inline void evaluate_residual(size_t level, std::vector<std::vector<T>> &result);

        /** @brief Compute L_{H}U_{H} on a given grid of resolution H */
        inline void evaluate_operator(size_t level, ndArray<T, Ndim> &result);
        /** @brief Compute compute the residual */
        inline void evaluate_residual(size_t level, ndArray<T, Ndim> &result);

        /** @brief Apply N times the smoother */
        void relax(const size_t level, size_t N);

        /** @brief Apply the smoother until convergence with respect to the given tolerance or until the maximum iteration */
        void relax(const size_t level, const double tolerance);

        /** @brief specialization class methods (e.g Liniar multigrid or NonLinear multigrid) */
        virtual inline void multigrid() {}
        virtual inline void solve() { multigrid(); }

        /** @brief get the differential operator at a given position in 2D*/  
        // specialized through class function specialization
        virtual T differential_operator(size_t, int, int) = 0;
        /** @brief get the differential operator at a given position in 3D*/
         // specialized through class function specialization
        virtual T differential_operator(size_t, int, int, int) = 0;
        /** @brief relation at a given position in 2D*/
         // specialized through class function specialization
        virtual void relaxation_updater(size_t, int, int) = 0;
        /** @brief relation at a given position in 3D*/
         // specialized through class function specialization
        virtual void relaxation_updater(size_t, int, int, int) = 0;

    protected:
        Stack<T, Ndim> solution, source;   // Grids for solution and source term
        Stack<T, Ndim> temp;               // Extra storage for muligrid solver
        CycleType cycle_type;              // cycle type (V, W or F)
        size_t Npre;                       // Number of pre-smoothing step
        size_t Npost;                      // Number of post-smoothing step
        const size_t maxIter;              // Maximum iteration
        const T residualTolerance;         // Tolerance
        size_t finestLevel, coarsestLevel; // finest and coarsest levels
        size_t nxfine, nyfine, nzfine;     // grid resolutions
        bool isSourceSet;                  // flag related to source setting
        bool isInitialValueSet;            // flag related to initial guess setting

    private:
        T residualSum, normSum; // residual and solution norms
        Derivatives du;
    };

    template <typename T, int Ndim>
    MultigridBase<T, Ndim>::MultigridBase(const Settings<T> &settings) :

                                                                         solution(settings),
                                                                         source(settings),
                                                                         temp(settings),
                                                                         cycle_type(settings.cycleType),
                                                                         Npre(settings.Npre),
                                                                         Npost(settings.Npost),
                                                                         residualTolerance(settings.residualTolerance),
                                                                         maxIter(settings.maxIterations),
                                                                         isSourceSet(false),
                                                                         isInitialValueSet(false)
    {
        finestLevel = solution.finestLevel;
        coarsestLevel = solution.coarsestLevel;
        if (Ndim == 2)
        {
            nxfine = solution[finestLevel].size(0);
            nyfine = solution[finestLevel].size(1);
        }
        else if (Ndim == 3)
        {
            nxfine = solution[finestLevel].size(0);
            nyfine = solution[finestLevel].size(1);
            nzfine = solution[finestLevel].size(2);
        }
    }

    // template <typename T, int Ndim>
    // inline std::vector<std::vector<T>> &MultigridBase<T, Ndim>::get_result()
    // {
    //     return solution[finestLevel];
    // }

    template <typename T, int Ndim>
    inline ndArray<T, Ndim> &MultigridBase<T, Ndim>::get_result()
    {
        return solution[finestLevel];
    }

    template <typename T, int Ndim>
    template <typename U>
    inline void MultigridBase<T, Ndim>::initial_guess(U arg)
    {
        solution[finestLevel] = arg;
        isInitialValueSet = true;
    }

    // template <typename T, int Ndim>
    // inline std::vector<std::vector<T>> &MultigridBase<T, Ndim>::source_term()
    // {
    //     return source[finestLevel];
    // }

    template <typename T, int Ndim>
    inline ndArray<T, Ndim> &MultigridBase<T, Ndim>::source_term()
    {
        return source[finestLevel];
    }

    template <typename T, int Ndim>
    template <typename U>
    inline void MultigridBase<T, Ndim>::source_term(U arg)
    {
        source[finestLevel] = arg;
        isSourceSet = true;
    }

    template <typename T, int Ndim>
    inline void MultigridBase<T, Ndim>::evaluate_operator(size_t level, std::vector<std::vector<T>> &result)
    {
        size_t nx = result.size();
        size_t ny = result[0].size();
        for (auto i = 0; i < nx; i++)
        {
            for (auto j = 0; j < ny; j++)
            {
                result[i][j] = differential_operator(level, i, j);
            }
        }
    }

    template <typename T, int Ndim>
    inline void MultigridBase<T, Ndim>::evaluate_residual(size_t level, std::vector<std::vector<T>> &result)
    {
        size_t nx = result.size();
        size_t ny = result[0].size();
        for (auto i = 0; i < nx; i++)
        {
            for (auto j = 0; j < ny; j++)
            {
                result[i][j] = source[level][i][j] - differential_operator(level, i, j);
            }
        }
    }

    template <typename T, int Ndim>
    inline void MultigridBase<T, Ndim>::evaluate_operator(size_t level, ndArray<T, Ndim> &result)
    {
        size_t nx = result.size(0);
        size_t ny = result.size(1);
        size_t nz = 0;
        if (Ndim == 3)
            nz = result.size(2);

        if (Ndim == 2)
        {
            for (auto i = 0; i < nx; i++)
            {
                for (auto j = 0; j < ny; j++)
                {
                    result({i, j}) = differential_operator(level, i, j);
                }
            }
        }

        else if (Ndim == 3)
        {
            for (auto k = 0; k < nz; k++)
            {
                for (auto j = 0; j < ny; j++)
                {
                    for (auto i = 0; i < nx; i++)
                    {
                        result({i, j, k}) = differential_operator(level, i, j, k);
                    }
                }
            }
        }
    }

    template <typename T, int Ndim>
    inline void MultigridBase<T, Ndim>::evaluate_residual(size_t level, ndArray<T, Ndim> &result)
    {
        size_t nx = result.size(0);
        size_t ny = result.size(1);
        size_t nz = 0;
        if (Ndim == 3)
            nz = result.size(2);

        if (Ndim == 2)
        {
            for (auto j = 0; j < ny; j++)
            {
                for (auto i = 0; i < nx; i++)
                {
                    result({i, j}) = source[level]({i, j}) - differential_operator(level, i, j);
                }
            }
        }

        else if (Ndim == 3)
        {
            for (auto k = 0; k < nz; k++)
            {
                for (auto j = 0; j < ny; j++)
                {
                    for (auto i = 0; i < nx; i++)
                    {
                        result({i, j, k}) = source[level]({i, j, k}) - differential_operator(level, i, j, k);
                    }
                }
            }
        }
    }


// RB-GS smoother
    template <typename T, int Ndim>
    void MultigridBase<T, Ndim>::relax(const size_t level, size_t N)
    {
        size_t nx{solution[level].size(0)}, ny{solution[level].size(1)};
        size_t nz = (Ndim == 3) ? solution[level].size(2) : 0;

        // Col-major order
        for (auto round = 0; round < N; round++)
        {
            for (size_t pColor = 2; pColor > 0; pColor--)
            {
                if (Ndim == 2)
                {
                    for (size_t j = 1; j < ny - 1; j++)
                    {
                        auto i = pColor - power(-1, pColor) * (j % 2);
                        for (; i < nx - 1; i += 2)
                        {
                            relaxation_updater(level, i, j);
                        }
                    }
                }

                else if (Ndim == 3)
                {
                    for (size_t k = 1; k < nz - 1; k++)
                    {
                        for (size_t j = 1; j < ny - 1; j++)
                        {
                            auto i = pColor - power(-1, pColor) * ((k + j) % 2);
                            for (; i < nx - 1; i += 2)
                            {
                                relaxation_updater(level, i, j, k);
                            }
                        }
                    }
                }
            }
            // TODO : Boundaries update
        }

        // Row-major order RB-GS
        /*
        for(size_t pColor = 2; pColor > 0; pColor--)
        {
            if(Ndim == 2)
            {
                for(size_t i = 1; i < nx - 1; i++)
                {
                    auto j = pColor + power(-1,pColor-1) * (i%2);
                    for(; j < ny -1; j+=2)
                    {
                        relaxation_updater(level, i, j);
                    }
                }
            }

            else if(Ndim == 3)
            {
                for(size_t i = 1; j < nx-1; i++)
                {
                    for(size_t j = 1; j < ny -1; j++)
                    {
                        auto k = pColor + power(-1, pColor) * ((i + j) % 2);
                        for(; k < nz - 1; k+=2)
                        {
                            relaxation_updater(level, i,j,k);
                        }
                    }
                }
            }
        }
        */
    }


// RB-GS smoother
    template <typename T, int Ndim>
    void MultigridBase<T, Ndim>::relax(const size_t level, const double tolerance)
    {
        size_t nx{solution[level].size(0)}, ny{solution[level].size(1)};
        size_t nz = (Ndim == 3) ? solution[level].size(2) : 0;

        for (size_t iter = 0; iter < maxIter; iter++)
        {
            T L2ResidualNorm = 0, L2SolutionNorm = 0, tmp_buf = 0;
            /*Col-major mode*/

            for (size_t pColor = 2; pColor > 0; pColor--)
            {
                if (Ndim == 2)
                {
                    for (size_t j = 1; j < ny - 1; j++)
                    {
                        auto i = pColor - power(-1, pColor) * (j % 2);
                        for (; i < nx - 1; i += 2)
                        {
                            tmp_buf = solution[level]({i, j});
                            relaxation_updater(level, i, j);

                            // current defect
                            // TODO include boundaries cells in the residuals norm
                            tmp_buf = solution[level]({i, j}) - tmp_buf;
                            L2ResidualNorm += power(tmp_buf, 2);
                            L2SolutionNorm += power(solution[level]({i, j}), 2);
                        }
                    }
                }

                else if (Ndim == 3)
                {
                    for (size_t k = 1; k < nz - 1; k++)
                    {
                        for (size_t j = 1; j < ny - 1; j++)
                        {
                            auto i = pColor - power(-1, pColor) * ((k + j) % 2);
                            for (; i < nx - 1; i += 2)
                            {
                                tmp_buf = solution[level]({i, j, k});
                                relaxation_updater(level, i, j, k);
                                // current defect
                                // TODO include boundaries cells in the residuals norm
                                tmp_buf = solution[level]({i, j, k}) - tmp_buf;
                                L2ResidualNorm += power(tmp_buf, 2);
                                L2SolutionNorm += power(solution[level]({i, j, k}), 2);
                            }
                        }
                    }
                }

                // TODO : Boundaries update

                if ((sqrt(L2ResidualNorm) / sqrt(L2SolutionNorm)) < tolerance)
                    return;
            }
        }
    }

}

#endif