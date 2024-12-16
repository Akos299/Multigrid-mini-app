#ifndef __MULTIGRID_HPP__
#define __MULTIGRID_HPP__

#include <cstddef>
#include <array>

#include "../include/utilities.hpp"
#include "../include/boundary_condition.hpp"
#include "../include/mgrid_stack.hpp"
#include "../ndarray/ndArray.hpp"

namespace multigrid
{

    template <typename T, int Ndim>
    T L2Norm(nd::ndArray<T, Ndim> &in_arr)
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
        MultigridBase(const Settings<T> &settings, std::array<double, 3> &starts, std::array<double, 3> &ends, TransfertOperator &rest_operator, TransfertOperator &inter_operator);
        virtual ~MultigridBase() {};

        /** @brief get the output of the multigrid method */
        inline nd::ndArray<T, Ndim> &get_solution();

        /** @brief set the initial guess (e.g potential in the case of Poisson equation) */
        template <typename U>
        inline void set_initial_guess(U arg);

        /** @brief get the soruce term */
        inline nd::ndArray<T, Ndim> &get_source();

        /** @brief set the source term with the given argument */
        template <typename U>
        inline void set_source(U arg);

        /** @brief Compute L_{H}U_{H} on a given grid of resolution H */
        inline void evaluate_operator(size_t level, nd::ndArray<T, Ndim> &result);

        /** @brief Compute compute the residual/defect */
        inline void evaluate_residual(size_t level, nd::ndArray<T, Ndim> &result);

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
        // virtual T differential_operator(size_t, int, int) = 0;

        /** @brief get the differential operator at a given position in 3D*/
        // specialized through class function specialization
        // virtual T differential_operator(size_t, int, int, int) = 0;

        /** @brief relation at a given position in 2D*/
        // specialized through class function specialization
        // virtual void relaxation_updater(size_t, int, int, T w = 1.0) = 0;
        void relaxation_updater(size_t, int, int, T w = 1.0);

        /** @brief relation at a given position in 3D*/
        // specialized through class function specialization
        // virtual void relaxation_updater(size_t, int, int, int, T w = 1.0) = 0;
        void relaxation_updater(size_t, int, int, int, T w = 1.0);

        /** @brief return the source at the current level*/
        nd::ndArray<T, Ndim> &get_current_source() { return source[current_level]; }

        /** @brief return the Solution (e.g potential) at the current level*/
        nd::ndArray<T, Ndim> &get_current_solution() { return solution[current_level]; }

        /** @brief return the old solution at the current*/
        nd::ndArray<T, Ndim> &get_current_old_solution() { return old_solution[current_level]; }

        /** @brief return the current level*/
        size_t get_current_level() { return current_level; }

        /** @brief return the total number of level */
        size_t get_total_number_of_level() { return nb_levels; }

        /** @brief restrict the solution throw all level
         *          This will be useful for the full multigrid.
         *        Endeed by restriction until the coarse level, one can then start the FMG algo.
         */
        inline void restrict_solution();

        /** @brief restrict the the source through all level*/
        inline void restrict_source();

        /** @brief set solution at the given level to zero*/
        inline void set_solution_to_zero(size_t level);

        /** @brief set source  at the given level to zero*/
        inline void set_source_to_zero(size_t level);

        /** @brief set solution at all level to zero  */
        inline void set_all_solution_to_zero();

        /** @brief set source at all level  to zero*/
        inline void set_all_source_to_zero();

        /** @brief perform One step down : Npre smooth + Restriction to the coarsest level*/
        inline void one_step_down();

        /** @brief perform One step up :  Prolongation + correction + Npost smooth*/
        inline void one_step_up();

        /** @brief perform V-cyle multigrid iteration starting from the current level */
        inline void mgi_Vcycle();

        /** @brief perform W-cycle multigrid iteration starting from the current level */
        inline void mgi_Wcycle();

        /** @brief perform F-cycle multigrid iteration starting from the current level*/
        inline void mgi_Fcycle();

        /** @brief perform multigrid iteration starting from the current level*/
        inline void mgi(const CycleType &cycle_type);

        /** @brief perform a cycle of full multigrid  */
        inline void fmg_cycle(const CycleType &cycle_type, const int mgc_per_level = 1.0);

        /** @brief perform a  of full multigrid  */
        inline void fmg(const CycleType &cycle_type, const int mgc_per_level = 1.);

        /** @brief perform iteration of V(or W or F)-cycle until convergence.
         * This can be done by specify a given threshold.
         */
        inline void iterate_mg_cycle_to_convergence(const CycleType &cycle_type, const T eps);

        /** @brief perform iteration of V(or W or F)-cycle until convergence.
         * This can be done by specify a number of iteration to saturate
         */
        inline void iterate_mg_cycle_to_convergence(const CycleType &cycle_type, const int n_iter);

    protected:
        Stack<T, Ndim> solution, source, old_solution; // Grids for solution and source term
        Stack<T, Ndim> temp;                           // Extra storage for muligrid solver
        CycleType cycle_type;                          // cycle type (V, W or F)
        size_t npre;                                   // Number of pre-smoothing step
        size_t npost;                                  // Number of post-smoothing step
        const size_t max_iter;                         // Maximum iteration
        const T res_tol;                               // Tolerance
        size_t finest_level, coarsest_level;           // finest and coarsest levels
        size_t nx_fine, ny_fine, nz_fine;              // grid (finest level) resolutions
        bool is_source_set;                            // flag related to source setting
        bool is_initial_value_set;                     // flag related to initial guess setting
        size_t current_level;                          // the current level
        size_t nb_levels;                              // Total number of level
        T w;                                           // smoothing parameter
        bool fixedNiter;
        std::array<double, 3> lim_start, lim_end;
        TransfertOperator res_ope, int_ope;
        BcType bc;

    private:
        T residual_sum, norm_sum; // residual and solution norms
        Derivatives du;
    };

}

#endif