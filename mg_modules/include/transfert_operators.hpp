#ifndef __TRANSFERT_OPERATORS__
#define __TRANSFERT_OPERATORS__

#include "../ndarray/ndArray.hpp"
#include "../include/utilities.hpp"

namespace multigrid
{

    /********************************************** Transfert opertors ************************************************
     *       interpolation operators  and restriction operators :
     *       To give the user the choice to experiment different transfert opertors we implement severals prolongation
     *       and restriction operators. Bellow are the list of the one available right now. But note that the list is eveoving.
     *                             @@@ Operators list  @@@
     *                              --------------------
     *         I - Prolongation : (bi/tri)-linear, Wesselling(or Khalil), Kwak, constantwise
     *         II - Restriction : (bi/tri)-linear, Wesselling(or Khalil), Kwak, constantwise
     *********************************************************************************************************************
     */
    /** @brief linear restriction operator */
    template <typename T, int Ndim>
    inline void linear_restriction_operator(nd::ndArray<T, Ndim> &corseData, nd::ndArray<T, Ndim> &fineData)
    {
        size_t nxc, nyc, nzc, nxf, nyf, nzf;
        nxc = corseData.size(0);
        nyc = corseData.size(1);
        nzc = (Ndim == 3) ? (corseData.size(2)) : 0;
        nxf = fineData.size(0);
        nyf = fineData.size(1);
        nzf = (Ndim == 3) ? (fineData.size(2)) : 0;

        if (Ndim == 2)
        {
            // update iner cells
            for (int j = 1; j < nyc - 1; j++)
            {
                auto fj = 2 * j;
                for (int i = 1; i < nxc - 1; i++)
                {
                    auto fi = 2 * i;
                    corseData({i, j}) = (1 / 64) * (9 * fineData({fi, fj}) + 9 * fineData({fi + 1, fj}) + 9 * fineData({fi + 1, fj + 1}) + 9 * fineData({fi, fj + 1}) + 3 * fineData({fi - 1, fj}) + 3 * fineData({fi + 2, fj}) + 3 * fineData({fi + 2, fj + 1}) + 3 * fineData({fi - 1, fj + 1}) + fineData({fi - 1, fj - 1}) + fineData({fi + 2, fj - 1}) + fineData({fi + 2, fj + 2}) + fineData({fi - 1, fj + 2}) + 3 * fineData({fi, fj - 1}) + 3 * fineData({fi + 1, 2 * j - 1}) + 3 * fineData({fi + 1, fj + 2}) + 3 * fineData({fi, fj + 2}));
                }
            }

            // update boundaries cells (deepending on the BC types e.g Dirichlet, Neumann, Robin, etc..)
        }

        else if (Ndim == 3)
        {
            // update iner cells
            for (int k = 1; k < nzc - 1; k++)
            {
                auto fk = 2 * k;
                for (int j = 1; j < nyc - 1; j++)
                {
                    auto fj = 2 * j;
                    for (int i = 1; i < nxc - 1; i++)
                    {
                        auto fi = 2 * i;
                        corseData({i, j, k}) = (1 / 512) * (fineData({fi - 1, fj + 2, fk - 1}) + 3 * fineData({fi - 1, fj + 1, fk - 1}) + 3 * fineData({fi - 1, fj, fk - 1}) + fineData({fi - 1, fj - 1, fk - 1}) + 3 * fineData({fi, fj + 2, fk - 1}) + 9 * fineData({fi, fj + 1, fk - 1}) + 9 * fineData({fi, fj, fk - 1}) + 3 * fineData({fi, fj - 1, fk - 1}) + 3 * fineData({fi + 1, fj + 2, fk - 1}) + 9 * fineData({fi + 1, fj + 1, fk - 1}) + 9 * fineData({fi + 1, fj, fk - 1}) + 3 * fineData({fi + 1, fj - 1, fk - 1}) + fineData({fi + 2, fj + 2, fk - 1}) + 3 * fineData({fi + 2, fj + 1, fk - 1}) + 3 * fineData({fi + 2, fj, fk - 1}) + fineData({fi + 2, fj - 1, fk - 1})

                                                            + 3 * fineData({fi - 1, fj + 2, fk}) + 9 * fineData({fi - 1, fj + 1, fk}) + 9 * fineData({fi - 1, fj, fk}) + 3 * fineData({fi - 1, fj - 1, fk}) + 9 * fineData({fi, fj + 2, fk}) + 27 * fineData({fi, fj + 1, fk}) + 27 * fineData({fi, fj, fk}) + 9 * fineData({fi, fj - 1, fk}) + 9 * fineData({fi + 1, fj + 2, fk}) + 27 * fineData({fi + 1, fj + 1, fk}) + 27 * fineData({fi + 1, fj, fk}) + 9 * fineData({fi + 1, fj - 1, fk}) + 3 * fineData({fi + 2, fj + 2, fk}) + 9 * fineData({fi + 2, fj + 1, fk}) + 9 * fineData({fi + 2, fj, fk}) + 3 * fineData({fi + 2, fj - 1, fk})

                                                            + 3 * fineData({fi - 1, fj + 2, fk + 1}) + 9 * fineData({fi - 1, fj + 1, fk + 1}) + 9 * fineData({fi - 1, fj, fk + 1}) + 3 * fineData({fi - 1, fj - 1, fk + 1}) + 9 * fineData({fi, fj + 2, fk + 1}) + 27 * fineData({fi, fj + 1, fk + 1}) + 27 * fineData({fi, fj, fk + 1}) + 9 * fineData({fi, fj - 1, fk + 1}) + 9 * fineData({fi + 1, fj + 2, fk + 1}) + 27 * fineData({fi + 1, fj + 1, fk + 1}) + 27 * fineData({fi + 1, fj, fk + 1}) + 9 * fineData({fi + 1, fj - 1, fk + 1}) + 3 * fineData({fi + 2, fj + 2, fk + 1}) + 9 * fineData({fi + 2, fj + 1, fk + 1}) + 9 * fineData({fi + 2, fj, fk + 1}) + 3 * fineData({fi + 2, fj - 1, fk + 1})

                                                            + fineData({fi - 1, fj + 2, fk + 2}) + 3 * fineData({fi - 1, fj + 1, fk + 2}) + 3 * fineData({fi - 1, fj, fk + 2}) + fineData({fi - 1, fj - 1, fk + 2}) + 3 * fineData({fi, fj + 2, fk + 2}) + 9 * fineData({fi, fj + 1, fk + 2}) + 9 * fineData({fi, fj, fk + 2}) + 3 * fineData({fi, fj - 1, fk + 2}) + 3 * fineData({fi + 1, fj + 2, fk + 2}) + 9 * fineData({fi + 1, fj + 1, fk + 2}) + 9 * fineData({fi + 1, fj, fk + 2}) + 3 * fineData({fi + 1, fj - 1, fk + 2}) + fineData({fi + 2, fj + 2, fk + 2}) + 3 * fineData({fi + 2, fj + 1, fk + 2}) + 3 * fineData({fi + 2, fj, fk + 2}) + fineData({fi + 2, fj - 1, fk + 2}));
                    }
                }
            }

            // updates boundaries cells
        }
    }

    /** @brief khalil restriction operator */
    template <typename T, int Ndim>
    inline void khalil_restriction_operator(nd::ndArray<T, Ndim> &corseData, nd::ndArray<T, Ndim> &fineData)
    {
        size_t nxc, nyc, nzc, nxf, nyf, nzf;
        nxc = corseData.size(0);
        nyc = corseData.size(1);
        nzc = (Ndim == 3) ? (corseData.size(2)) : 0;
        nxf = fineData.size(0);
        nyf = fineData.size(1);
        nzf = (Ndim == 3) ? (fineData.size(2)) : 0;

        if (Ndim == 2)
        {
            // update iner cells
            for (int j = 1; j < nyc - 1; j++)
            {
                auto fj = 2 * j;
                for (int i = 1; i < nxc - 1; i++)
                {
                    auto fi = 2 * i;
                    corseData({i, j}) = (1 / 16) * (2 * fineData({fi, fj}) + 3 * fineData({fi + 1, fj}) + 3 * fineData({fi, fj + 1}) + 3 * fineData({fi + 1, fj + 1}) + fineData({fi, fj + 2}) + 0 * fineData({fi, fj - 1}) + 0 * fineData({fi + 1, fj + 2}) + fineData({fi + 1, fj - 1}) + 0 * fineData({fi - 1, fj}) + fineData({fi - 1, fj + 1}) + 0 * fineData({fi + 2, fj + 1}) + 1 * fineData({fi + 2, fj}) + fineData({fi - 1, fj + 2}) + 0 * fineData({fi - 1, fj - 1}) + 0 * fineData({fi + 2, fj + 2}) + fineData({fi + 1, fj - 1}));
                }
            }

            // update boundaries cells (deepending on the BC types e.g Dirichlet, Neumann, Robin, etc..)
        }

        else if (Ndim == 3)
        {
            // update iner cells
            for (int k = 1; k < nzc - 1; k++)
            {
                auto fk = 2 * k;
                for (int j = 1; j < nyc - 1; j++)
                {
                    auto fj = 2 * j;
                    for (int i = 1; i < nxc - 1; i++)
                    {
                        auto fi = 2 * i;
                        corseData({i, j, k}) = (1 / 32) * (fineData({fi + 1, fj, fk - 1}) + fineData({fi + 1, fj - 1, fk - 1}) + fineData({fi + 2, fj, fk - 1}) + fineData({fi + 2, fj - 1, fk - 1}) + 2 * fineData({fi, fj + 1, fk}) + 2 * fineData({fi, fj, fk}) + 2 * fineData({fi + 1, fj + 1, fk}) + 3 * fineData({fi + 1, fj, fk}) + fineData({fi + 1, fj - 1, fk}) + fineData({fi + 2, fj, fk}) + fineData({fi + 2, fj - 1, fk}) + fineData({fi - 1, fj + 2, fk + 1}) + fineData({fi - 1, fj + 1, fk + 1}) + fineData({fi, fj + 2, fk + 1}) + 3 * fineData({fi, fj + 1, fk + 1}) + 2 * fineData({fi, fj, fk + 1}) + 2 * fineData({fi + 1, fj + 1, fk + 1}) + 2 * fineData({fi + 1, fj, fk + 1}) + fineData({fi - 1, fj + 2, fk + 2}) + fineData({fi - 1, fj + 1, fk + 2}) + fineData({fi, fj + 2, fk + 2}) + fineData({fi, fj + 1, fk + 2}));
                    }
                }
            }

            // updates boundaries cells
        }
    }

    /** @brief kwak restriction operator */
    template <typename T, int Ndim>
    inline void kwak_restriction_operator(nd::ndArray<T, Ndim> &corseData, nd::ndArray<T, Ndim> &fineData)
    {
        size_t nxc, nyc, nzc, nxf, nyf, nzf;
        nxc = corseData.size(0);
        nyc = corseData.size(1);
        nzc = (Ndim == 3) ? (corseData.size(2)) : 0;
        nxf = fineData.size(0);
        nyf = fineData.size(1);
        nzf = (Ndim == 3) ? (fineData.size(2)) : 0;

        if (Ndim == 2)
        {
            // update iner cells
            for (int j = 1; j < nyc - 1; j++)
            {
                auto fj = 2 * j;
                for (int i = 1; i < nxc - 1; i++)
                {
                    auto fi = 2 * i;
                    corseData({i, j}) = (1 / 16) * (2 * fineData({fi, fj}) + 2 * fineData({fi + 1, fj}) + 2 * fineData({fi, fj + 1}) + 2 * fineData({fi + 1, fj + 1}) + fineData({fi, fj + 2}) + fineData({fi, fj - 1}) + fineData({fi + 1, fj + 2}) + fineData({fi + 1, fj - 1}) + fineData({fi - 1, fj}) + fineData({fi - 1, fj + 1}) + fineData({fi + 2, fj + 1}) + fineData({fi + 2, fj}) + 0 * fineData({fi - 1, fj + 2}) + 0 * fineData({fi - 1, fj - 1}) + 0 * fineData({fi + 2, fj + 2}) + 0 * fineData({fi + 1, fj - 1}));
                }
            }

            // update boundaries cells (deepending on the BC types e.g Dirichlet, Neumann, Robin, etc..)
        }

        else if (Ndim == 3)
        {
            // update iner cells
            for (int k = 1; k < nzc - 1; k++)
            {
                auto fk = 2 * k;
                for (int j = 1; j < nyc - 1; j++)
                {
                    auto fj = 2 * j;
                    for (int i = 1; i < nxc - 1; i++)
                    {
                        auto fi = 2 * i;
                        corseData({i, j, k}) = (1 / 48) * (fineData({i, j + 1, k - 1}) + fineData({i, j, k - 1}) + fineData({i + 1, j, k - 1}) + fineData({i + 1, j + 1, k - 1})

                                                           + fineData({i - 1, j + 1, k}) + fineData({i - 1, j, k}) + fineData({i, j + 2, k}) + 3 * fineData({i, j + 1, k}) + 3 * fineData({i, j, k}) + fineData({i, j - 1, k}) + fineData({i + 1, j + 2, k}) + 3 * fineData({i + 1, j + 1, k}) + 3 * fineData({i + 1, j, k}) + fineData({i + 1, j - 1, k}) + fineData({i + 2, j + 1, k}) + fineData({i + 2, j, k})

                                                           + fineData({i - 1, j + 1, k + 1}) + fineData({i - 1, j, k + 1}) + fineData({i, j + 2, k + 1}) + 3 * fineData({i, j + 1, k + 1}) + 3 * fineData({i, j, k + 1}) + fineData({i, j - 1, k + 1}) + fineData({i + 1, j + 2, k + 1}) + 3 * fineData({i + 1, j + 1, k + 1}) + 3 * fineData({i + 1, j, k + 1}) + fineData({i + 1, j - 1, k + 1}) + fineData({i + 2, j + 1, k + 1}) + fineData({i + 2, j, k + 1})

                                                           + fineData({i, j + 1, k + 2}) + fineData({i, j, k + 2}) + fineData({i + 1, j, k + 2}) + fineData({i + 1, j + 1, k + 2}));
                    }
                }
            }

            // updates boundaries cells
        }
    }

    /** @brief constant restriction operator */
    template <typename T, int Ndim>
    inline void constantwise_restriction_operator(nd::ndArray<T, Ndim> &corseData, nd::ndArray<T, Ndim> &fineData)
    {
        size_t nxc, nyc, nzc, nxf, nyf, nzf;
        nxc = corseData.size(0);
        nyc = corseData.size(1);
        nzc = (Ndim == 3) ? (corseData.size(2)) : 0;
        nxf = fineData.size(0);
        nyf = fineData.size(1);
        nzf = (Ndim == 3) ? (fineData.size(2)) : 0;

        if (Ndim == 2)
        {
            // update iner cells
            for (int j = 1; j < nyc - 1; j++)
            {
                auto fj = 2 * j;
                for (int i = 1; i < nxc - 1; i++)
                {
                    auto fi = 2 * i;
                    corseData({i, j}) = (1 / 4) * (fineData({fi, fj}) + fineData({fi + 1, fj}) + fineData({fi, fj + 1}) + fineData({fi + 1, fj + 1}));
                }
            }

            // update boundaries cells (deepending on the BC types e.g Dirichlet, Neumann, Robin, etc..)
        }

        else if (Ndim == 3)
        {
            // update iner cells
            for (int k = 1; k < nzc - 1; k++)
            {
                auto fk = 2 * k;
                for (int j = 1; j < nyc - 1; j++)
                {
                    auto fj = 2 * j;
                    for (int i = 1; i < nxc - 1; i++)
                    {
                        auto fi = 2 * i;
                        corseData({i, j, k}) = (1 / 8) * (fineData({fi, fj, fk}) + fineData({fi + 1, fj, fk}) + fineData({fi, fj + 1, fk}) + fineData({fi, fj, fk + 1}) + fineData({fi + 1, fj, fk + 1}) + fineData({fi, fj + 1, fk + 1}) + fineData({fi + 1, fj + 1, fk}) + fineData({fi + 1, fj + 1, fk + 1}));
                    }
                }
            }

            // updates boundaries cells
        }
    }

    /** @brief constantwise prolongation operator */
    template <typename T, int Ndim>
    inline void constantwise_prolongation_operator(nd::ndArray<T, Ndim> &corseData, nd::ndArray<T, Ndim> &fineData)
    {
        size_t nxc, nyc, nzc, nxf, nyf, nzf;
        nxc = corseData.size(0);
        nyc = corseData.size(1);
        nzc = (Ndim == 3) ? (corseData.size(2)) : 0;
        nxf = fineData.size(0);
        nyf = fineData.size(1);
        nzf = (Ndim == 3) ? (fineData.size(2)) : 0;

        if (Ndim == 2)
        {
            // update iner cells
            for (int j = 0; j < nyc; j++)
            {
                auto fj = 2 * j;
                for (int i = 0; i < nxc; i++)
                {
                    auto fi = 2 * i;
                    auto tmp = corseData({i, j});
                    fineData({fi, fj}) = tmp;
                    fineData({fi + 1, fj}) = tmp;
                    fineData({fi + 1, fj + 1}) = tmp;
                    fineData({fi, fj + 1}) = tmp;
                }
            }
        }

        if (Ndim == 3)
        {
            // update iner cells
            for (int k = 0; k < nzc; k++)
            {
                auto fk = 2 * k;
                for (int j = 0; j < nyc; j++)
                {
                    auto fj = 2 * j;
                    for (int i = 0; i < nxc; i++)
                    {
                        auto fi = 2 * i;
                        auto tmp = corseData({i, j, k});
                        fineData({fi, fj, fk}) = tmp;
                        fineData({fi + 1, fj, fk}) = tmp;
                        fineData({fi + 1, fj + 1, fk}) = tmp;
                        fineData({fi, fj + 1, fk}) = tmp;

                        fineData({fi, fj, fk + 1}) = tmp;
                        fineData({fi + 1, fj, fk + 1}) = tmp;
                        fineData({fi + 1, fj + 1, fk + 1}) = tmp;
                        fineData({fi, fj + 1, fk + 1}) = tmp;
                    }
                }
            }
        }
    }

    /** @brief linear prolongation operator */
    // TODO make sure nd::ndArray is initialize to zero
    template <typename T, int Ndim>
    inline void linear_prolongation_operator(nd::ndArray<T, Ndim> &corseData, nd::ndArray<T, Ndim> &fineData)
    {
        size_t nxc, nyc, nzc, nxf, nyf, nzf;
        nxc = corseData.size(0);
        nyc = corseData.size(1);
        nzc = (Ndim == 3) ? (corseData.size(2)) : 0;
        nxf = fineData.size(0);
        nyf = fineData.size(1);
        nzf = (Ndim == 3) ? (fineData.size(2)) : 0;

        if (Ndim == 2)
        {
            // update iner cells
            for (int j = 1; j < nyc - 1; j++)
            {
                auto fj = 2 * j;
                for (int i = 1; i < nxc - 1; i++)
                {
                    auto fi = 2 * i;
                    fineData({fi, fj}) = (1 / 16) * (9 * corseData({i, j}) + 3 * corseData({i - 1, j}) + 3 * corseData({i, j - 1}) + corseData({i - 1, j - 1}));
                    fineData({fi + 1, fj}) = (1 / 16) * (9 * corseData({i, j}) + 3 * corseData({i + 1, j}) + 3 * corseData({i, j - 1}) + corseData({i + 1, j - 1}));
                    fineData({fi, fj + 1}) = (1 / 16) * (9 * corseData({i, j}) + 3 * corseData({i, j + 1}) + 3 * corseData({i - 1, j}) + corseData({i - 1, j + 1}));
                    fineData({fi + 1, fj + 1}) = (1 / 16) * (9 * corseData({i, j}) + 3 * corseData({i, j + 1}) + 3 * corseData({i + 1, j}) + corseData({i + 1, j + 1}));
                }
            }

            // TODO BC to add on
        }
        if (Ndim == 3)
        {
            for (int k = 1; k < nzc - 1; k++)
            {

                auto fk = 2 * k;
                for (int j = 1; j < nyc - 1; j++)
                {
                    auto fj = 2 * j;
                    for (int i = 1; i < nxc - 1; i++)
                    {
                        auto fi = 2 * i;

                        fineData({fi, fj, fk}) = (1 / 64) * (27 * corseData({i, j, k}) + 9 * corseData({i - 1, j, k}) + 3 * corseData({i - 1, j - 1, k}) + 9 * corseData({i, j - 1, k}) + 3 * corseData({i - 1, j, k - 1}) + corseData({i - 1, j - 1, k - 1}) + 3 * corseData({i, j - 1, k - 1}) + 9 * corseData({i, j, k - 1}));

                        fineData({fi + 1, fj, fk}) = (1 / 64) * (27 * corseData(i, j, k) + 9 * corseData({i, j - 1, k}) + 3 * corseData({i + 1, j - 1, k}) + 9 * corseData(i + 1, j, k) + 3 * corseData({i, j - 1, k - 1}) + corseData({i + 1, j - 1, k - 1}) + 9 * corseData({i, j, k - 1}) + 3 * corseData({i + 1, j, k - 1}));

                        fineData({fi, fj + 1, fk}) = (1 / 64) * (27 * corseData(i, j, k) + 9 * corseData({i, j + 1, k}) + 9 * corseData({i - 1, j, k}) + 3 * corseData({i - 1, j + 1, k}) + 9 * corseData({i, j, k - 1}) + 3 * corseData({i, j + 1, k - 1}) + 3 * corseData({i - 1, j, k - 1}) + corseData({i - 1, j + 1, k - 1}));

                        fineData({fi + 1, fj + 1, fk}) = (1 / 64) * (27 * corseData(i, j, k) + 9 * corseData({i + 1, j, k}) + 9 * corseData({i, j + 1, k}) + 3 * corseData(i + 1, j + 1, k) + 9 * corseData({i, j, k - 1}) + 3 * corseData({i + 1, j, k - 1}) + 3 * corseData({i, j + 1, k - 1}) + corseData(i + 1, j + 1, k - 1));

                        fineData({fi, fj, fk + 1}) = (1 / 64) * (27 * corseData({i, j, k}) + 9 * corseData({i, j - 1, k}) + 9 * corseData({i - 1, j, k}) + 3 * corseData({i - 1, j - 1, k}) + 9 * corseData({i, j, k + 1}) + 3 * corseData({i, j - 1, k + 1}) + 3 * corseData({i - 1, j, k + 1}) + corseData({i - 1, j - 1, k + 1}));

                        fineData({fi + 1, fj, fk + 1}) = (1 / 64) * (27 * corseData({i, j, k}) + 9 * corseData({i + 1, j, k}) + 9 * corseData({i, j - 1, k}) + 3 * corseData(i + 1, j - 1, k) + 9 * corseData({i, j, k + 1}) + 3 * corseData({i + 1, j, k + 1}) + 3 * corseData(i, j - 1, k + 1) + corseData({i + 1, j - 1, k + 1}));

                        fineData({fi, fj + 1, fk + 1}) = (1 / 64) * (27 * corseData({i, j, k}) + 9 * corseData({i, j + 1, k}) + 9 * corseData(i - 1, j, k) + 3 * corseData(i - 1, j + 1, k) + 9 * corseData({i, j, k + 1}) + 3 * corseData({i - 1, j, k + 1}) + 3 * corseData({i, j + 1, k + 1}) + corseData({i - 1, j + 1, k + 1}));

                        fineData({fi + 1, fj + 1, fk + 1}) = (1 / 64) * (27 * corseData({i, j, k}) + 9 * corseData({i + 1, j, k}) + 9 * corseData({i, j + 1, k}) + corseData({i + 1, j + 1, k}) + 9 * corseData({i, j, k + 1}) + 3 * corseData({i, j + 1, k + 1}) + 3 * corseData({i + 1, j, k + 1}) + corseData({i + 1, j + 1, k + 1}));
                    }
                }
            }

            // TODO : add BC
        }
    }

    /** @brief Khalil prolongation operator */
    // TODO make sure nd::ndArray is initialize to zero
    template <typename T, int Ndim>
    inline void khalil_prolongation_operator(nd::ndArray<T, Ndim> &corseData, nd::ndArray<T, Ndim> &fineData)
    {
        size_t nxc, nyc, nzc, nxf, nyf, nzf;
        nxc = corseData.size(0);
        nyc = corseData.size(1);
        nzc = (Ndim == 3) ? (corseData.size(2)) : 0;
        nxf = fineData.size(0);
        nyf = fineData.size(1);
        nzf = (Ndim == 3) ? (fineData.size(2)) : 0;

        if (Ndim == 2)
        {
            // update iner cells
            for (int j = 1; j < nyc - 1; j++)
            {
                auto fj = j * 2;

                for (int i = 1; i < nxc - 1; i++)
                {
                    auto fi = 2 * i;
                    fineData({fi, fj}) = (1 / 4) * (2 * corseData({i, j}));
                    fineData({fi + 1, fj}) = (1 / 4) * (3 * corseData({i, j}) + corseData({i + 1, j}) + corseData({i, j - 1}) + corseData({i + 1, j - 1}));
                    fineData({fi, fj + 1}) = (1 / 4) * (3 * corseData({i, j}) + corseData({i - 1, j}) + corseData({i, j + 1}) + corseData({i - 1, j + 1}));
                    fineData({fi + 1, fj + 1}) = (1 / 4) * (2 * corseData(i, j));
                }
            }

            // TODO BC to add on
        }

        if (Ndim == 3)
        {
            for (int k = 1; k < nzc - 1; k++)
            {
                auto fk = 2 * k;
                for (int j = 1; j < nyc - 1; j++)
                {
                    auto fj = 2 * j;
                    for (int i = 1; i < nxc - 1; i++)
                    {
                        auto fi = 2 * i;

                        fineData({fi, fj, fk}) = (1 / 4) * (2 * corseData({i, j, k}));

                        fineData({fi + 1, fj, fk}) = (1 / 4) * (3 * corseData(i, j, k) + corseData({i, j - 1, k}) + corseData({i + 1, j - 1, k}) + corseData(i + 1, j, k) + corseData({i, j - 1, k - 1}) + corseData({i + 1, j - 1, k - 1}) + corseData({i, j, k - 1}) + corseData({i + 1, j, k - 1}));

                        fineData({fi, fj + 1, fk}) = (1 / 4) * (2 * corseData(i, j, k));

                        fineData({fi + 1, fj + 1, fk}) = (1 / 4) * (2 * corseData(i, j, k));

                        fineData({fi, fj, fk + 1}) = (1 / 4) * (2 * corseData({i, j, k}));

                        fineData({fi + 1, fj, fk + 1}) = (1 / 4) * (2 * corseData({i, j, k}));

                        fineData({fi, fj + 1, fk + 1}) = (1 / 4) * (3 * corseData({i, j, k}) + corseData({i, j + 1, k}) + corseData(i - 1, j, k) + corseData(i - 1, j + 1, k) + corseData({i, j, k + 1}) + corseData({i - 1, j, k + 1}) + corseData({i, j + 1, k + 1}) + corseData({i - 1, j + 1, k + 1}));

                        fineData({fi + 1, fj + 1, fk + 1}) = (1 / 4) * (2 * corseData({i, j, k}));
                    }
                }
            }

            // TODO : add BC
        }
    }

    /** @brief Kwak prolongation operator */
    // TODO make sure nd::ndArray is initialize to zero
    template <typename T, int Ndim>
    inline void kwak_prolongation_operator(nd::ndArray<T, Ndim> &corseData, nd::ndArray<T, Ndim> &fineData)
    {
        size_t nxc, nyc, nzc, nxf, nyf, nzf;
        nxc = corseData.size(0);
        nyc = corseData.size(1);
        nzc = (Ndim == 3) ? (corseData.size(2)) : 0;
        nxf = fineData.size(0);
        nyf = fineData.size(1);
        nzf = (Ndim == 3) ? (fineData.size(2)) : 0;

        if (Ndim == 2)
        {
            // update iner cells
            for (int j = 1; j < nyc - 1; j++)
            {
                auto fj = 2 * j;
                for (int i = 1; i < nxc - 1; i++)
                {
                    auto fi = 2 * i;
                    fineData({fi, fj}) = (1 / 4) * (2 * corseData({i, j}) + corseData({i - 1, j}) + corseData({i, j - 1}));
                    fineData({fi + 1, fj}) = (1 / 4) * (2 * corseData({i, j}) + corseData({i, j - 1}) + corseData({i + 1, j}));
                    fineData({fi, fj + 1}) = (1 / 4) * (2 * corseData({i, j}) + corseData({i - 1, j}) + corseData({i, j + 1}));
                    fineData({fi + 1, fj + 1}) = (1 / 4) * (2 * corseData({i, j}) + corseData({i + 1, j}) + corseData({i + 1, j + 1}));
                }
            }

            // TODO BC to add on
        }

        if (Ndim == 3)
        {
            for (int k = 1; k < nzc - 1; k++)
            {
                auto fk = 2 * k;
                for (int j = 1; j < nyc - 1; j++)
                {
                    auto fj = 2 * j;
                    for (int i = 1; i < nxc - 1; i++)
                    {
                        auto fi = 2 * i;

                        fineData({fi, fj, fk}) = (1 / 6) * (3 * corseData({i, j, k}) + corseData({i - 1, j, k}) + corseData({i - 1, j - 1, k}) + corseData({i, j, k - 1}));

                        fineData({fi + 1, fj, fk}) = (1 / 6) * (3 * corseData({i, j, k}) + corseData({i, j - 1, k}) + corseData({i + 1, j, k}) + corseData({i, j, k - 1}));

                        fineData({fi, fj + 1, fk}) = (1 / 6) * (3 * corseData({i, j, k}) + corseData({i, j + 1, k}) + corseData({i - 1, j, k}) + corseData({i, j, k - 1}));

                        fineData({fi + 1, fj + 1, fk}) = (1 / 6) * (3 * corseData(i, j, k) + corseData({i + 1, j, k}) + corseData({i, j + 1, k}) + corseData({i, j, k - 1}));

                        fineData({fi, fj, fk + 1}) = (1 / 6) * (3 * corseData({i, j, k}) + corseData({i, j - 1, k}) + corseData({i - 1, j, k}) + corseData({i, j, k + 1}));

                        fineData({fi + 1, fj, fk + 1}) = (1 / 6) * (3 * corseData({i, j, k}) + corseData({i + 1, j, k}) + corseData({i, j - 1, k}) + corseData({i, j, k + 1}));

                        fineData({fi, fj + 1, fk + 1}) = (1 / 6) * (3 * corseData({i, j, k}) + corseData({i, j + 1, k}) + corseData(i - 1, j, k) + corseData({i, j, k + 1}));

                        fineData({fi + 1, fj + 1, fk + 1}) = (1 / 6) * (3 * corseData({i, j, k}) + corseData({i + 1, j, k}) + corseData({i, j + 1, k}) + corseData({i, j, k + 1}));
                    }
                }
            }

            // update BC
        }
    }

    template <typename T, TransfertOperator Op, int Ndim>
    inline void cubic_operator(nd::ndArray<T, Ndim> &corseData, nd::ndArray<T, Ndim> &fineData)
    {
        size_t nxc, nyc, nzc, nxf, nyf, nzf;
        nxc = corseData.size(0);
        nyc = corseData.size(1);
        nzc = (Ndim == 3) ? (corseData.size(2)) : 0;
        nxf = fineData.size(0);
        nyf = fineData.size(1);
        nzf = (Ndim == 3) ? (fineData.size(2)) : 0;

        if (Ndim == 2)
        {
            std::cout << "cubic interpolation not available in 2D. Use bilinear interpolator." << "\n";
        }

        if (Ndim == 3)
        {

            for (int k = 1; k < nzc - 1; k++)
            {

                auto fk = 2 * k;
                for (int j = 1; j < nyc - 1; j++)
                {
                    auto fj = 2 * j;
                    for (int i = 1; i < nxc - 1; i++)
                    {
                        auto fi = 2 * i;

                        fineData({fi, fj, fk}) = (1.0 / 32768.0) * (125. * corseData({i - 1, j - 1, k - 1}) + 750. * corseData({i, j - 1, k - 1}) - 75. * corseData({i + 1, j - 1, k - 1}) + 750. * corseData({i - 1, j, k - 1}) + 4500. * corseData({i, j, k - 1}) - 450. * corseData({i + 1, j, k - 1}) - 75. * corseData({i - 1, j + 1, k - 1}) - 450. * corseData({i, j + 1, k - 1}) + 45. * corseData({i + 1}, j + 1, k - 1) + 750. * corseData({i - 1, j - 1, k}) + 4500. * corseData({i, j - 1, k}) - 450. * corseData({i + 1, j - 1, k}) + 4500. * corseData({i - 1, j, k}) + 27000. * corseData({i, j, k}) - 2700. * corseData({i + 1, j, k}) - 450. * corseData({i - 1, j + 1, k}) - 2700. * corseData({i, j + 1, k}) + 270. * corseData({i + 1, j + 1, k}) - 75. * corseData({i - 1, j - 1, k + 1}) - 450. * corseData({i, j - 1, k + 1}) + 45. * corseData({i + 1, j - 1, k + 1}) - 450. * corseData({i - 1, j, k + 1}) - 2700. * corseData({i, j, k + 1}) + 270. * corseData({i + 1, j, k + 1}) + 45. * corseData({i - 1, j + 1, k + 1}) + 270. * corseData({i, j + 1, k + 1}) - 27. * corseData({i + 1, j + 1, k + 1}));

                        fineData({fi + 1, fj, fk}) = (1.0 / 32768.0) * (-75. * corseData({i - 1, j - 1, k - 1}) + 750. * corseData({i, j - 1, k - 1}) + 125. * corseData({i + 1, j - 1, k - 1}) - 450. * corseData({i - 1, j, k - 1}) + 4500. * corseData({i, j, k - 1}) + 750. * corseData({i + 1, j, k - 1}) + 45. * corseData({i - 1, j + 1, k - 1}) - 450. * corseData({i, j + 1, k - 1}) - 75. * corseData({i + 1}, j + 1, k - 1) - 450. * corseData({i - 1, j - 1, k}) + 4500. * corseData({i, j - 1, k}) + 750. * corseData({i + 1, j - 1, k}) - 2700. * corseData({i - 1, j, k}) + 27000. * corseData({i, j, k}) + 4500. * corseData({i + 1, j, k}) + 270. * corseData({i - 1, j + 1, k}) - 2700. * corseData({i, j + 1, k}) - 450. * corseData({i + 1, j + 1, k}) + 45. * corseData({i - 1, j - 1, k + 1}) - 450. * corseData({i, j - 1, k + 1}) - 75. * corseData({i + 1, j - 1, k + 1}) + 270. * corseData({i - 1, j, k + 1}) - 2700. * corseData({i, j, k + 1}) - 450. * corseData({i + 1, j, k + 1}) - 27. * corseData({i - 1, j + 1, k + 1}) + 270. * corseData({i, j + 1, k + 1}) + 45. * corseData({i + 1, j + 1, k + 1}));

                        fineData({fi, fj + 1, fk}) = (1.0 / 32768.0) * (-75. * corseData({i - 1, j - 1, k - 1}) - 450. * corseData({i, j - 1, k - 1}) + 45. * corseData({i + 1, j - 1, k - 1}) + 750. * corseData({i - 1, j, k - 1}) + 4500. * corseData({i, j, k - 1}) - 450. * corseData({i + 1, j, k - 1}) + 125. * corseData({i - 1, j + 1, k - 1}) + 750. * corseData({i, j + 1, k - 1}) - 75. * corseData({i + 1}, j + 1, k - 1) - 450. * corseData({i - 1, j - 1, k}) - 2700. * corseData({i, j - 1, k}) + 270. * corseData({i + 1, j - 1, k}) + 4500. * corseData({i - 1, j, k}) + 27000. * corseData({i, j, k}) - 2700. * corseData({i + 1, j, k}) + 750. * corseData({i - 1, j + 1, k}) + 4500. * corseData({i, j + 1, k}) - 450. * corseData({i + 1, j + 1, k}) + 45. * corseData({i - 1, j - 1, k + 1}) + 270. * corseData({i, j - 1, k + 1}) - 27. * corseData({i + 1, j - 1, k + 1}) - 450. * corseData({i - 1, j, k + 1}) - 2700. * corseData({i, j, k + 1}) + 270. * corseData({i + 1, j, k + 1}) - 75. * corseData({i - 1, j + 1, k + 1}) - 450. * corseData({i, j + 1, k + 1}) + 45. * corseData({i + 1, j + 1, k + 1}));

                        fineData({fi + 1, fj + 1, fk}) = (1.0 / 32768.0) * (+45. * corseData({i - 1, j - 1, k - 1}) - 450. * corseData({i, j - 1, k - 1}) - 75. * corseData({i + 1, j - 1, k - 1}) - 450. * corseData({i - 1, j, k - 1}) + 4500. * corseData({i, j, k - 1}) + 750. * corseData({i + 1, j, k - 1}) - 75. * corseData({i - 1, j + 1, k - 1}) + 750. * corseData({i, j + 1, k - 1}) + 125. * corseData({i + 1}, j + 1, k - 1) + 270. * corseData({i - 1, j - 1, k}) - 2700. * corseData({i, j - 1, k}) - 450. * corseData({i + 1, j - 1, k}) - 2700. * corseData({i - 1, j, k}) + 27000. * corseData({i, j, k}) + 4500. * corseData({i + 1, j, k}) - 450. * corseData({i - 1, j + 1, k}) + 4500. * corseData({i, j + 1, k}) + 750. * corseData({i + 1, j + 1, k}) - 27. * corseData({i - 1, j - 1, k + 1}) + 270. * corseData({i, j - 1, k + 1}) + 45. * corseData({i + 1, j - 1, k + 1}) + 270. * corseData({i - 1, j, k + 1}) - 2700. * corseData({i, j, k + 1}) - 450. * corseData({i + 1, j, k + 1}) + 45. * corseData({i - 1, j + 1, k + 1}) - 450. * corseData({i, j + 1, k + 1}) - 75. * corseData({i + 1, j + 1, k + 1}));

                        fineData({fi, fj, fk + 1}) = (1.0 / 32768.0) * (-75. * corseData({i - 1, j - 1, k - 1}) - 450. * corseData({i, j - 1, k - 1}) + 45. * corseData({i + 1, j - 1, k - 1}) - 450. * corseData({i - 1, j, k - 1}) - 2700. * corseData({i, j, k - 1}) + 270. * corseData({i + 1, j, k - 1}) + 45. * corseData({i - 1, j + 1, k - 1}) + 270. * corseData({i, j + 1, k - 1}) - 27. * corseData({i + 1}, j + 1, k - 1) + 750. * corseData({i - 1, j - 1, k}) + 4500. * corseData({i, j - 1, k}) - 450. * corseData({i + 1, j - 1, k}) + 4500. * corseData({i - 1, j, k}) + 27000. * corseData({i, j, k}) - 2700. * corseData({i + 1, j, k}) - 450. * corseData({i - 1, j + 1, k}) - 2700. * corseData({i, j + 1, k}) + 270. * corseData({i + 1, j + 1, k}) + 125. * corseData({i - 1, j - 1, k + 1}) + 750. * corseData({i, j - 1, k + 1}) - 75. * corseData({i + 1, j - 1, k + 1}) + 750. * corseData({i - 1, j, k + 1}) + 4500. * corseData({i, j, k + 1}) - 450. * corseData({i + 1, j, k + 1}) - 75. * corseData({i - 1, j + 1, k + 1}) - 450. * corseData({i, j + 1, k + 1}) + 45. * corseData({i + 1, j + 1, k + 1}));

                        fineData({fi + 1, fj, fk + 1}) = (1.0 / 32768.0) * (+45. * corseData({i - 1, j - 1, k - 1}) - 450. * corseData({i, j - 1, k - 1}) - 75. * corseData({i + 1, j - 1, k - 1}) + 270. * corseData({i - 1, j, k - 1}) - 2700. * corseData({i, j, k - 1}) - 450. * corseData({i + 1, j, k - 1}) - 27. * corseData({i - 1, j + 1, k - 1}) + 270. * corseData({i, j + 1, k - 1}) + 45. * corseData({i + 1}, j + 1, k - 1) - 450. * corseData({i - 1, j - 1, k}) + 4500. * corseData({i, j - 1, k}) + 750. * corseData({i + 1, j - 1, k}) - 2700. * corseData({i - 1, j, k}) + 27000. * corseData({i, j, k}) + 4500. * corseData({i + 1, j, k}) + 270. * corseData({i - 1, j + 1, k}) - 2700. * corseData({i, j + 1, k}) - 450. * corseData({i + 1, j + 1, k}) - 75. * corseData({i - 1, j - 1, k + 1}) + 750. * corseData({i, j - 1, k + 1}) + 125. * corseData({i + 1, j - 1, k + 1}) - 450. * corseData({i - 1, j, k + 1}) + 4500. * corseData({i, j, k + 1}) + 750. * corseData({i + 1, j, k + 1}) + 45. * corseData({i - 1, j + 1, k + 1}) - 450. * corseData({i, j + 1, k + 1}) - 75. * corseData({i + 1, j + 1, k + 1}));

                        fineData({fi, fj + 1, fk + 1}) = (1.0 / 32768.0) * (+45. * corseData({i - 1, j - 1, k - 1}) + 270. * corseData({i, j - 1, k - 1}) - 27. * corseData({i + 1, j - 1, k - 1}) - 450. * corseData({i - 1, j, k - 1}) - 2700. * corseData({i, j, k - 1}) + 270. * corseData({i + 1, j, k - 1}) - 75. * corseData({i - 1, j + 1, k - 1}) - 450. * corseData({i, j + 1, k - 1}) + 45. * corseData({i + 1}, j + 1, k - 1) - 450. * corseData({i - 1, j - 1, k}) - 2700. * corseData({i, j - 1, k}) + 270. * corseData({i + 1, j - 1, k}) + 4500. * corseData({i - 1, j, k}) + 27000. * corseData({i, j, k}) - 2700. * corseData({i + 1, j, k}) + 750. * corseData({i - 1, j + 1, k}) + 4500. * corseData({i, j + 1, k}) - 450. * corseData({i + 1, j + 1, k}) - 75. * corseData({i - 1, j - 1, k + 1}) - 450. * corseData({i, j - 1, k + 1}) + 45. * corseData({i + 1, j - 1, k + 1}) + 750. * corseData({i - 1, j, k + 1}) + 4500. * corseData({i, j, k + 1}) - 450. * corseData({i + 1, j, k + 1}) + 125. * corseData({i - 1, j + 1, k + 1}) + 750. * corseData({i, j + 1, k + 1}) - 75. * corseData({i + 1, j + 1, k + 1}));

                        fineData({fi + 1, fj + 1, fk + 1}) = (1.0 / 32768.0) * (-27. * corseData({i - 1, j - 1, k - 1}) + 270. * corseData({i, j - 1, k - 1}) + 45. * corseData({i + 1, j - 1, k - 1}) + 270. * corseData({i - 1, j, k - 1}) - 2700. * corseData({i, j, k - 1}) + 450. * corseData({i + 1, j, k - 1}) + 45. * corseData({i - 1, j + 1, k - 1}) - 450. * corseData({i, j + 1, k - 1}) - 75. * corseData({i + 1}, j + 1, k - 1) + 270. * corseData({i - 1, j - 1, k}) - 2700. * corseData({i, j - 1, k}) - 450. * corseData({i + 1, j - 1, k}) - 2700. * corseData({i - 1, j, k}) + 27000. * corseData({i, j, k}) + 4500. * corseData({i + 1, j, k}) + 450. * corseData({i - 1, j + 1, k}) + 4500. * corseData({i, j + 1, k}) + 750. * corseData({i + 1, j + 1, k}) + 45. * corseData({i - 1, j - 1, k + 1}) - 450. * corseData({i, j - 1, k + 1}) - 75. * corseData({i + 1, j - 1, k + 1}) - 450. * corseData({i - 1, j, k + 1}) + 4500. * corseData({i, j, k + 1}) + 750. * corseData({i + 1, j, k + 1}) - 75. * corseData({i - 1, j + 1, k + 1}) + 750. * corseData({i, j + 1, k + 1}) + 125. * corseData({i + 1, j + 1, k + 1}));
                    }
                }
            }

            // TODO : add BC
        }
    }

    /** @brief Te;plate based prolongation operator
     */
    template <typename T, TransfertOperator Op, int Ndim>
    inline void prolongation_operator(nd::ndArray<T, Ndim> &corseData, nd::ndArray<T, Ndim> &fineData)
    {
        if constexpr (Op == PointAverage)
            constantwise_prolongation_operator<T, Ndim>(corseData, fineData);
        if constexpr (Op == Linear)
            linear_prolongation_operator<T, Ndim>(corseData, fineData);
        if constexpr (Op == Khalil)
            khalil_prolongation_operator<T, Ndim>(corseData, fineData);
        if constexpr (Op == Kwak)
            kwak_prolongation_operator<T, Ndim>(corseData, fineData);
        if constexpr (Op == TreeCubic)
            cubic_operator<T, Ndim>(corseData, fineData);
    }

    /** @brief Template based restriction operator
     */
    template <typename T, TransfertOperator Op, int Ndim>
    inline void restriction_operator(nd::ndArray<T, Ndim> &corseData, nd::ndArray<T, Ndim> &fineData)
    {
        if constexpr (Op == PointAverage)
            constantwise_restriction_operator<T, Ndim>(corseData, fineData);
        if constexpr (Op == Linear)
            linear_restriction_operator<T, Ndim>(corseData, fineData);
        if constexpr (Op == Khalil)
            khalil_restriction_operator<T, Ndim>(corseData, fineData);
        if constexpr (Op == Kwak)
            kwak_restriction_operator(corseData, fineData);
    }
}

#endif