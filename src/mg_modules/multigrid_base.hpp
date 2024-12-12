#ifndef __MULTIGRID__
#define __MULTIGRID__

#include <cstddef>
#include <vector>
#include <array>
#include "utilities.hpp"
#include "boundary_condition.hpp"
#include "mgrid_stack.hpp"
#include "ndArray.h"

namespace multigrid
{

    template <typename T, int Ndim>
    T L2Norm(ndArray<T, Ndim> &in_arr)
    {
        T res = 0.;
        size_t nx = in_arr.size(0), ny = in_arr.size(1), nz = in_arr.size(2);
        if (Ndim == 2)
        {
            for (int j = 0; j < ny; j++)
            {
                for (int i = 0; i < nx; i++)
                {
                    res += power(in_arr({i, j}), 2);
                }
            }
        }
        if (Ndim == 3)
        {
            for (int k = 0; k < nz; k++)
            {
                for (int j = 0; j < ny; j++)
                {
                    for (int i = 0; i < nx; i++)
                    {
                        res += power(in_arr({i, j, k}), 2);
                    }
                }
            }
        }
        return sqrt(res);
    }

    /*********************** MultigridBase class******************
     *  This is a class that's will be specialized for each type of problem
     */
    template <typename T, int Ndim>
    class MultigridBase
    {
    public:
        MultigridBase(const Settings<T> &settings);
        virtual ~MultigridBase() {};

        /** @brief get the output of the multigrid method */
        inline ndArray<T, Ndim> &getSolution();

        /** @brief set the initial guess (e.g potential in the case of Poisson equation) */
        template <typename U>
        inline void setInitialGuess(U arg);

        /** @brief get the soruce term */
        inline ndArray<T, Ndim> &getSource();

        /** @brief set the source term with the given argument */
        template <typename U>
        inline void setSource(U arg);

        /** @brief Compute L_{H}U_{H} on a given grid of resolution H */
        inline void evaluateOperator(size_t level, ndArray<T, Ndim> &result);

        /** @brief Compute compute the residual/defect */
        inline void evaluateResidual(size_t level, ndArray<T, Ndim> &result);

        /** @brief Apply N times the smoother */
        template <DiscreteOperator Lh>
        void smoother(const size_t level, size_t N, T w = 1.0);

        /** @brief Apply the smoother until convergence with respect to the given tolerance or until the maximum iteration
         *           Use this to get the coarse level solution
         */
        template <DiscreteOperator Lh>
        void smoother(const size_t level, const double tolerance, T w = 1.0);

        /** @brief specialization class methods (e.g Liniar multigrid or NonLinear multigrid) */
        virtual inline void multigrid() {}

        virtual inline void solve() { multigrid(); }

        /** @brief get the differential operator at a given position in 2D*/
        // specialized through class function specialization
        virtual T differentialOperator(size_t, int, int) = 0;

        /** @brief get the differential operator at a given position in 3D*/
        // specialized through class function specialization
        virtual T differentialOperator(size_t, int, int, int) = 0;

        /** @brief relation at a given position in 2D*/
        // specialized through class function specialization
        virtual void relaxationUpdater(size_t, int, int, T w = 1.0) = 0;

        /** @brief relation at a given position in 3D*/
        // specialized through class function specialization
        virtual void relaxationUpdater(size_t, int, int, int, T w = 1.0) = 0;

        /** @brief return the source at the current level*/
        ndArray<T, Ndim> &getCurrentSource() { return source[current_level]; }

        /** @brief return the Solution (e.g potential) at the current level*/
        ndArray<T, Ndim> &getCurrentSolution() { return solution[current_level]; }

        /** @brief return the old solution at the current*/
        ndArray<T, Ndim> &getCurrentOldSolution() { return oldSolution[current_level]; }

        /** @brief return the current level*/
        size_t getCurrentLevel() { return current_level; }

        /** @brief return the total number of level */
        size_t getTotalNumberOfLevel() { return nLevel; }

        /** @brief restrict the solution throw all level
         *          This will be useful for the full multigrid.
         *        Endeed by restriction until the coarse level, one can then start the FMG algo.
         */
        inline void restrictSolution();

        /** @brief restrict the the source through all level*/
        inline void restrictSource();

        /** @brief set solution at the given level to zero*/
        inline void setSolutioToZero(size_t level);

        /** @brief set source  at the given level to zero*/
        inline void setSourceToZero(size_t level);

        /** @brief set solution at all level to zero  */
        inline void setAllSolutionToZero();

        /** @brief set source at all level  to zero*/
        inline void setAllSourceToZero();

        /** @brief perform One step down : Npre smooth + Restriction to the coarsest level*/
        inline void oneStepDown();

        /** @brief perform One step up :  Prolongation + correction + Npost smooth*/
        inline void oneStepUp();

        /** @brief perform V-cyle multigrid iteration starting from the current level */
        inline void mgiVcycle();

        /** @brief perform W-cycle multigrid iteration starting from the current level */
        inline void mgiWcycle();

        /** @brief perform F-cycle multigrid iteration starting from the current level*/
        inline void mgiFcycle();

        /** @brief perform multigrid iteration starting from the current level*/
        inline void mgi(const CycleType &cycle_type);

        /** @brief perform a cycle of full multigrid  */
        inline void fmgCycle(const CycleType &cycle_type, const int mgc_per_level = 1.0);

        /** @brief perform a  of full multigrid  */
        inline void fmg(const CycleType &cycle_type, const int mgc_per_level = 1.);

        /** @brief perform iteration of V(or W or F)-cycle until convergence.
         * This can be done by specify a given threshold.
         */
        inline void iterateMGCycleToConvergence(const CycleType &cycle_type, const T eps);

        /** @brief perform iteration of V(or W or F)-cycle until convergence.
         * This can be done by specify a number of iteration to saturate
         */
        inline void iterateMGCycleToConvergence(const CycleType &cycle_type, const int n_iter);

    protected:
        Stack<T, Ndim> solution, source, old_solution; // Grids for solution and source term
        Stack<T, Ndim> temp;                          // Extra storage for muligrid solver
        CycleType cycle_type;                         // cycle type (V, W or F)
        size_t npre;                                  // Number of pre-smoothing step
        size_t npost;                                 // Number of post-smoothing step
        const size_t max_iter;                         // Maximum iteration
        const T res_tol;                    // Tolerance
        size_t finest_level, coarsest_level;            // finest and coarsest levels
        size_t nx_fine, ny_fine, nz_fine;                // grid (finest level) resolutions
        bool is_source_set;                             // flag related to source setting
        bool is_initial_value_set;                       // flag related to initial guess setting
        size_t current_level;                         // the current level
        size_t nb_levels;                                // Total number of level
        T w;                                          // smoothing parameter
        bool fixedNiter;

    private:
        T residual_sum, norm_sum; // residual and solution norms
        Derivatives du;
    };

    template <typename T, int Ndim>
    MultigridBase<T, Ndim>::MultigridBase(const Settings<T> &settings) :

                                                                         solution(settings),
                                                                         source(settings),
                                                                         temp(settings),
                                                                         cycle_type(settings.cycle_type),
                                                                         Npre(settings.Npre),
                                                                         Npost(settings.Npost),
                                                                         residualTolerance(settings.residualTolerance),
                                                                         maxIter(settings.maxIterations),
                                                                         isSourceSet(false),
                                                                         isInitialValueSet(false),
                                                                         w(settings.w),
                                                                         fixedNiter(settings.fixedNiter)
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
    inline ndArray<T, Ndim> &MultigridBase<T, Ndim>::getSolution()
    {
        return solution[finestLevel];
    }

    template <typename T, int Ndim>
    template <typename U>
    inline void MultigridBase<T, Ndim>::setInitialGuess(U arg)
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
    inline ndArray<T, Ndim> &MultigridBase<T, Ndim>::getSource()
    {
        return source[finestLevel];
    }

    template <typename T, int Ndim>
    template <typename U>
    inline void MultigridBase<T, Ndim>::setSource(U arg)
    {
        source[finestLevel] = arg;
        isSourceSet = true;
    }

    template <typename T, int Ndim>
    inline void MultigridBase<T, Ndim>::evaluateOperator(size_t level, ndArray<T, Ndim> &result)
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
                    result({i, j}) = differentialOperator(level, i, j);
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
                        result({i, j, k}) = differentialOperator(level, i, j, k);
                    }
                }
            }
        }
    }

    template <typename T, int Ndim>
    inline void MultigridBase<T, Ndim>::evaluateResidual(size_t level, ndArray<T, Ndim> &result)
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
                    result({i, j}) = source[level]({i, j}) - differentialOperator(level, i, j);
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
                        result({i, j, k}) = source[level]({i, j, k}) - differentialOperator(level, i, j, k);
                    }
                }
            }
        }
    }

    // RB-GS smoother  TODO extended it to add JACOBI, SOR, RKL
    // Change the interface so it can take as argument an operator Lh <depending of the problem both in (2D/3D)>,
    // the unknown vector (e.g defect or potential) and the  RHS vector (e.g the source)
    // The name should be change to smoother

    template <typename T, int Ndim>
    template <DiscreteOperator Lh>
    void MultigridBase<T, Ndim>::smoother(const size_t level, size_t N, T w)
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
                            relaxationUpdater(level, i, j, w);
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
                                relaxationUpdater(level, i, j, k, w);
                            }
                        }
                    }
                }
            }
            // TODO : Boundaries update
        }
    }

    // RB-GS smoother
    template <typename T, int Ndim>
    template <DiscreteOperator Lh>
    void MultigridBase<T, Ndim>::smoother(const size_t level, const double tolerance, T w)
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
                            relaxationUpdater(level, i, j, w);

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
                                relaxationUpdater(level, i, j, k, w);
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

    template <typename T, int Ndim>
    inline void MultigridBase<T, Ndim>::restrictSolution()
    {
        for (auto level = finestLevel; level > coarsestLevel /*we prefer to use coarsestLevel insteed of 0*/; level--)
        {
            solution.coarsen(level);
        }
    }

    template <typename T, int Ndim>
    inline void MultigridBase<T, Ndim>::restrictSource()
    {
        for (auto level = finestLevel; level > coarsestLevel /*we prefer to use coarsestLevel insteed of 0*/; level--)
        {
            source.coarsen(level);
        }
    }

    template <typename T, int Ndim>
    inline void MultigridBase<T, Ndim>::setSolutioToZero(size_t level) { solution[level].set_zero(); }

    template <typename T, int Ndim>
    inline void MultigridBase<T, Ndim>::setAllSolutionToZero()
    {
        for (auto level = finestLevel; level >= 0; level--)
            setSolutioToZero(level);
    }

    template <typename T, int Ndim>
    inline void MultigridBase<T, Ndim>::setSourceToZero(size_t level) { source[level].set_zero(); }

    template <typename T, int Ndim>
    inline void MultigridBase<T, Ndim>::setAllSourceToZero()
    {
        for (auto level = finestLevel; level >= 0; level--)
            setSourceToZero(level);
    }

    template <typename T, int Ndim>
    inline void MultigridBase<T, Ndim>::oneStepDown()
    {
        if (current_level == coarsestLevel)
            return;

        // perform Npre smoothing step
        smoother(current_level, Npre, w);

        // evaluate the defect/residual in the current stage
        evaluateResidual(current_level, temp[current_level]);

        // restrict the defect to the next level. The restrict defect is store in the next level source ndArray, so it can be used as source for the next level step
        temp[current_level].coarsen(current_level, source[current_level - 1]);

        // set the initial guess for the next level to zero
        setSolutioToZero(current_level - 1);
        current_level--;
    }

    template <typename T, int Ndim>
    inline void MultigridBase<T, Ndim>::oneStepUp()
    {
        // prolongate the correction to the finest level
        solution[current_level].refine(current_level, temp[current_level + 1]);

        // apply correction
        solution[current_level + 1] += temp[current_level + 1];

        // perform Npost relaxation
        smoother(current_level + 1, Npost, w);
        current_level++;
    }

    template <typename T, int Ndim>
    inline void MultigridBase<T, Ndim>::mgiVcycle()
    {
        size_t startLevel = current_level;
        // Downstroke
        while (current_level != coarsestLevel)
        {
            oneStepDown();
        }
        // solve the coarsest level
        smoother(current_level, residualTolerance, w);
        // Upstroke
        while (current_level != startLevel)
        {
            oneStepUp();
        }
    }

    template <typename T, int Ndim>
    inline void MultigridBase<T, Ndim>::mgi(const CycleType &cycle_type)
    {
        if (cycle_type == CycleType::fCycle)
        {
            mgiFcycle();
        }
        else if (cycle_type == CycleType::vCycle)
        {
            mgiVcycle();
        }
        else if (cycle_type == CycleType::wCycle)
        {
            mgiWcycle();
        }
        else
        {
            std::cout << "user-defined multigrid cycle not implemented " << std::endl;
        }
    }

    template <typename T, int Ndim>
    inline void MultigridBase<T, Ndim>::fmgCycle(const CycleType &cycle_type, const int mgc_per_level = 1.0 /*Number of multigrid cycle per level*/)
    {
        for (int level = 0; level < finestLevel; level++)
        {
            solution.refine<cubic_operator>(level);
            current_level = level + 1;
            for (int i = 0; i < mgc_per_level; i++)
                mgi<T, Ndim>(cycle_type);
        }
    }

    template <typename T, int Ndim>
    inline void MultigridBase<T, Ndim>::fmg(const CycleType &cycle_type, const int mgc_per_level = 1.0 /*Number of multigrid cycle per level*/)
    {
        // perform a restriction from the finest level to the coarsest one
        restrictSolution<T, Ndim>();
        restrictSource<T, Ndim>();
        // solve the coarsest level
        smoother<T, Ndim>(coarsestLevel, residualTolerance, w);
        // perform the full multigrid cycle
        fmgCycle<T, Ndim>(cycle_type, mgc_per_level);

        /**
         * If in the full multigrid algorithm we apply one V(or W or F)-cycle
         * per level, then at the end of the full multigrid algorithm we need to
         * apply additional V(or W or F)-cycle iterations to obtain a fully-converged solution.
         *
         * See also remark 2.6.2 from Trottenberg et al.2001
         */
        if (fixedNiter)
            iterateMGCycleToConvergence<T, Ndim>(cycle_type, residualTolerance);
        else
            iterateMGCycleToConvergence<T, Ndim>(cycle_type, maxIter);
    }

    template <typename T, int Ndim>
    inline void MultigridBase<T, Ndim>::iterateMGCycleToConvergence(const CycleType &cycle_type, const T eps)
    {

        ndArray<T, Ndim> _result; /** size of the grid to be specified*/
        int n_iter = 0;

        // computed the residual or the defect ====> should be modified to implement the equation (20) int Tomida & Stone (2023)
        evaluateResidual<T, Ndim>(finestLevel, _result);

        // compute the norm of the residual or the defect

        T def_norm = L2Norm<T, Ndim>(_result);
        while (def_norm > eps)
        {
            mgi<T, Ndim>(cycle_type);
            evaluateResidual<T, Ndim>(finestLevel, _result);
            T def_norm_new = L2Norm<T, Ndim>(_result);

            if (def_norm_new / def_norm > 1.0)
            {
                std::cout << "MultigridBase<T,Ndim>::iterateMGCycleToConvergence is diverging" << std::endl;
                break;
            }
            if (n_iter > 100)
            {
                std::cout << "MultigridBase<T,Ndim>::iterateMGCycleToConvergence  take too much time : " << std::endl;
                std::cout << "Niter done :  " << n_iter << std::endl;
                break;
            }
            n_iter++;
        }
    }

    template <typename T, int Ndim>
    inline void MultigridBase<T, Ndim>::iterateMGCycleToConvergence(const CycleType &cycle_type, const int n_iter)
    {

        ndArray<T, Ndim> _result; /** size of the grid to be specified*/

        for (int i = 0; i < n_iter; i++)
        {
            mgi<T, Ndim>(cycle_type);
            // computed the residual or the defect ====> should be modified to implement the equation (20) int Tomida & Stone (2023)
            evaluateResidual<T, Ndim>(finestLevel, _result);
            T def_norm_new = L2Norm<T, Ndim>(_result);

            // TODO add in the function signature a std::vector<T> to store the history of the residual
        }
    }

}

#endif