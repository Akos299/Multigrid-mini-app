#ifndef __MULTIGRID_HPP__
#define __MULTIGRID_HPP__

#include <array>
#include <cstddef>

#include "../include/boundary_condition.hpp"
#include "../include/mgrid_stack.hpp"
#include "../include/ndArray.h"
#include "../include/utilities.hpp"

namespace multigrid {

template <class T, int Ndim> inline T L2Norm(nd::ndArray<T, Ndim> &in_arr) {
  T res = 0.;
  size_t nx = in_arr.size(0), ny = in_arr.size(1), nz = in_arr.size(2);
  if (Ndim == 2) {
    for (nd::index_t j = 0; j < ny; j++) {
      for (nd::index_t i = 0; i < nx; i++) {
        res += power(in_arr({i, j}), 2);
      }
    }
  }
  if (Ndim == 3) {
    for (nd::index_t k = 0; k < nz; k++) {
      for (nd::index_t j = 0; j < ny; j++) {
        for (nd::index_t i = 0; i < nx; i++) {
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
template <class T, int Ndim> class MultigridBase {
public:
  // MultigridBase(Settings<T> &settings);
  MultigridBase(my_settings &settings, std::array<double, 3> &starts,
                std::array<double, 3> &ends, TransfertOperator &rest_operator,
                TransfertOperator &inter_operator)
      : solution(settings), source(settings), temp(settings),
        old_solution(settings), lim_start(starts), lim_end(ends),
        res_ope(rest_operator), int_ope(inter_operator),
        cycle_type(settings.cycle_type), npre(settings.npre),
        npost(settings.npost), res_tol(settings.res_tol),
        max_iter(settings.max_iter), is_source_set(false),
        is_initial_value_set(false), w(settings.w),
        fixedNiter(settings.fixedNiter), nb_levels(settings.nb_levels),
        bc(settings.bc_condition_type) {

    std::cout << "I'm in the MultigridBase constructor" << "\n";
    finest_level = solution.get_finer_level();
    coarsest_level = solution.get_coarser_level();

    std::cout << "fin_lev = " << finest_level << " coa_lev = " << coarsest_level
              << "\n";
    if (Ndim == 2) {
      nx_fine = solution[finest_level].size(0);
      ny_fine = solution[finest_level].size(1);
    } else if (Ndim == 3) {
      nx_fine = solution[finest_level].size(0);
      ny_fine = solution[finest_level].size(1);
      nz_fine = solution[finest_level].size(2);
    }
    std::cout << " nx_fine = " << nx_fine << " ny_fine = " << ny_fine
              << "nz_fine = " << nz_fine << "\n";

    //=========  set Physical boundary conditions on the finest level==========
    if (bc == BcType::Periodic) {
      set_periodic_bc<T, Ndim>(solution[finest_level]);
    }

    else if (bc == BcType::ZeroFixed) {
      set_zero_fixed_bc(solution[finest_level]);
    }

    else if (bc == BcType::ZeroGradient) {
      set_zero_gradient_bc(solution[finest_level]);

    } else if (bc == BcType::Isolated) {

    }

    else {
      std::cout << "## User defined boundary conditions " << std::endl;
    }

    //  ======= set Boundary conditions on the coarser level =========
    for (int level = finest_level - 1; level >= 0 /* coarsest_level*/; level--)
      set_boundary_condition<T, Ndim>(solution[level],
                                      bc); // we used fixed-zero condition for
                                           // intermediate level for the moment
  };

  virtual ~MultigridBase() {};

  friend class Stack<T, Ndim>;

  /** @brief get the output of the multigrid method */
  inline nd::ndArray<T, Ndim> &get_solution();

  /** @brief set the initial guess (e.g potential in the case of Poisson
  equation) */
  template <typename U> inline void set_initial_guess(U arg);

  /** @brief get the soruce term */
  inline nd::ndArray<T, Ndim> &get_source();

  /** @brief set the source term with the given argument */
  template <typename U> inline void set_source(U arg);

  /** @brief Compute L_{H}U_{H} on a given grid of resolution H */
  inline void evaluate_operator(size_t level, nd::ndArray<T, Ndim> &result);

  /** @brief Compute compute the residual/defect */
  inline void evaluate_residual(size_t level, nd::ndArray<T, Ndim> &result);

  /** @brief Apply N times the smoother */
  //   template <DiscreteOperator Lh>
  void smoother(const size_t level, size_t N, T w = 1.0);

  /** @brief Apply the smoother until convergence with respect to the given
  tolerance or until the maximum iteration
   *           Use this to get the coarse level solution
   */
  //   template <DiscreteOperator Lh>
  void smoother(const size_t level, const double tolerance, T w = 1.0);

  // /** @brief specialization class methods (e.g Liniar multigrid or NonLinear
  // multigrid) */ virtual inline void multigrid() {}

  // virtual inline void solve() { multigrid(); }

  // /** @brief get the differential operator at a given position in 2D*/
  // // specialized through class function specialization
  // // virtual T differential_operator(size_t, int, int) = 0;

  // /** @brief get the differential operator at a given position in 3D*/
  // // specialized through class function specialization
  // // virtual T differential_operator(size_t, int, int, int) = 0;

  /** @brief relation at a given position in 2D*/
  // specialized through class function specialization
  // virtual void relaxation_updater(size_t, int, int, T w = 1.0) = 0;
  void relaxation_updater(size_t, nd::index_t, nd::index_t, T w = 1.0);

  /** @brief relation at a given position in 3D*/
  // specialized through class function specialization
  // virtual void relaxation_updater(size_t, int, int, int, T w = 1.0) = 0;
  void relaxation_updater(size_t, nd::index_t, nd::index_t, nd::index_t,
                          T w = 1.0);

  /** @brief return the source at the current level*/
  nd::ndArray<T, Ndim> &get_current_source() { return source[current_level]; }

  /** @brief return the Solution (e.g potential) at the current level*/
  nd::ndArray<T, Ndim> &get_current_solution() {
    return solution[current_level];
  }

  /** @brief return the old solution at the current*/
  nd::ndArray<T, Ndim> &get_current_old_solution() {
    return old_solution[current_level];
  }

  /** @brief return the current level*/
  size_t get_current_level() { return current_level; }

  /** @brief return the total number of level */
  size_t get_total_number_of_level() { return nb_levels; }

  /** @brief restrict the solution throw all level
   *          This will be useful for the full multigrid.
   *        Endeed by restriction until the coarse level, one can then start
   the FMG algo.
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

  /** @brief perform One step down : Npre smooth + Restriction to the
  coarsest level*/
  inline void one_step_down();

  /** @brief perform One step up :  Prolongation + correction + Npost
  smooth*/
  inline void one_step_up();

  /** @brief perform V-cyle multigrid iteration starting from the current
  level */
  inline void mgi_Vcycle();

  /** @brief perform W-cycle multigrid iteration starting from the current
  level */
  inline void mgi_Wcycle();

  /** @brief perform F-cycle multigrid iteration starting from the current
  level*/
  inline void mgi_Fcycle();

  /** @brief perform multigrid iteration starting from the current level*/
  inline void mgi(const CycleType &cycle_type);

  /** @brief perform a cycle of full multigrid  */
  inline void fmg_cycle(const CycleType &cycle_type,
                        const int mgc_per_level = 1.0);

  /** @brief perform a  of full multigrid  */
  inline void fmg(const CycleType &cycle_type, const int mgc_per_level = 1.);

  /** @brief perform iteration of V(or W or F)-cycle until convergence.
   * This can be done by specify a given threshold.
   */
  inline void iterate_mg_cycle_to_convergence_tol(const CycleType &cycle_type,
                                                  const T eps);

  /** @brief perform iteration of V(or W or F)-cycle until convergence.
   * This can be done by specify a number of iteration to saturate
   */
  inline void iterate_mg_cycle_to_convergence_niter(const CycleType &cycle_type,
                                                    const int n_iter);

protected:
  Stack<T, Ndim> solution, source,
      old_solution;                    // Grids for solution and source term
  Stack<T, Ndim> temp;                 // Extra storage for muligrid solver
  CycleType cycle_type;                // cycle type (V, W or F)
  size_t npre;                         // Number of pre-smoothing step
  size_t npost;                        //    Number of post-smoothing step
  const size_t max_iter;               // Maximum  iteration
  const T res_tol;                     // Tolerance
  size_t finest_level, coarsest_level; // finest and coarsest levels
  size_t nx_fine, ny_fine, nz_fine;    // grid (finest level) resolutions
  bool is_source_set;                  // flag related to
                                       /* source setting */
  bool is_initial_value_set;           // flag related to initial guess setting
  size_t current_level;                // the current level
  size_t nb_levels;                    // Total number of level
  T w;                                 // smoothing parameter
  bool fixedNiter;
  std::array<double, 3> lim_start, lim_end;
  TransfertOperator res_ope, int_ope;
  BcType bc;

  // private:
  //     T residual_sum, norm_sum; // residual and solution norms
  //     Derivatives du;
};

/**/
template <typename T, int Ndim>
void MultigridBase<T, Ndim>::relaxation_updater(size_t level, nd::index_t i,
                                                nd::index_t j, T w) {
  size_t nx{solution[level].size(0)}, ny{solution[level].size(1)}, nz{1};

  double hx, hy, hz, ihxx, ihyy, ihzz;

  hx = (lim_end[0] - lim_start[0]) / nx;
  hy = (lim_end[1] - lim_start[1]) / ny;

  ihxx = 1.0 / (hx * hx);
  ihyy = 1.0 / (hy * hy);

  auto sol_tmp =
      ((solution[level]({i - 1, j}) + solution[level]({i + 1, j})) * ihxx +
       (solution[level]({i, j - 1}) + solution[level]({i, j + 1})) * ihyy -
       source[level]({i, j})) /
      (2 * (ihxx + ihyy));
  solution[level]({i, j}) += w * sol_tmp;
}

/**/
template <typename T, int Ndim>
void MultigridBase<T, Ndim>::relaxation_updater(size_t level, nd::index_t i,
                                                nd::index_t j, nd::index_t k,
                                                T w) {
  size_t nx{solution[level].size(0)}, ny{solution[level].size(1)}, nz{1};

  double hx, hy, hz, ihxx, ihyy, ihzz;

  hx = (lim_end[0] - lim_start[0]) / nx;
  hy = (lim_end[1] - lim_start[1]) / ny;
  hz = (lim_end[2] - lim_start[2]) / nz;
  ihxx = 1.0 / (hx * hx);
  ihyy = 1.0 / (hy * hy);
  ihzz = 1.0 / (hz * hz);
  auto sol_tmp =
      ((solution[level]({i - 1, j, k}) + solution[level]({i + 1, j, k})) *
           ihxx +
       (solution[level]({i, j - 1, k}) + solution[level]({i, j + 1, k})) *
           ihyy +
       (solution[level]({i, j, k - 1}) + solution[level]({i, j, k + 1})) *
           ihzz -
       source[level]({i, j, k})) /
      (2 * (ihxx + ihyy + ihzz));
  solution[level]({i, j, k}) += w * sol_tmp;
}

/**/
template <typename T, int Ndim>
inline nd::ndArray<T, Ndim> &MultigridBase<T, Ndim>::get_solution() {
  return solution[finest_level];
}

/**/
template <typename T, int Ndim>
template <typename U>
inline void MultigridBase<T, Ndim>::set_initial_guess(U arg) {
  solution[finest_level] = arg;
  is_initial_value_set = true;
}

/**/
template <typename T, int Ndim>
inline nd::ndArray<T, Ndim> &MultigridBase<T, Ndim>::get_source() {
  return source[finest_level];
}

/**/
template <typename T, int Ndim>
template <typename U>
inline void MultigridBase<T, Ndim>::set_source(U arg) {
  source[finest_level] = arg;
  is_source_set = true;
}

/**/
template <typename T, int Ndim>
inline void
MultigridBase<T, Ndim>::evaluate_operator(size_t level,
                                          nd::ndArray<T, Ndim> &result) {
  size_t nx = result.size(0);
  size_t ny = result.size(1);
  size_t nz = 1;
  size_t i_start{1}, i_end{nx - 1}, j_start{1}, j_end{ny - 1}, k_start, k_end;
  double hx, hy, hz, ihxx, ihyy, ihzz;

  if (Ndim == 3) {
    nz = result.size(2);
    k_start = 1;
    k_end = nz - 1;
  }

  hx = (lim_end[0] - lim_start[0]) / nx;
  hy = (lim_end[1] - lim_start[1]) / ny;
  hz = (lim_end[2] - lim_start[2]) / nz;
  ihxx = 1.0 / (hx * hx);
  ihyy = 1.0 / (hy * hy);
  ihzz = 1.0 / (hz * hz);

  if (Ndim == 2) {
    for (auto i = i_start; i < i_end; i++) {
      for (auto j = j_start; j < j_end; j++) {
        // result({i, j}) = differential_operator(level, i, j);
        result({i, j}) = ihxx * (solution[level]({i - 1, j}) -
                                 2.0 * solution[level]({i, j}) +
                                 solution[level]({i + 1, j})) +
                         ihyy * (solution[level]({i, j - 1}) -
                                 2.0 * solution[level]({i, j}) +
                                 solution[level]({i, j + 1}));
      }
    }
  }

  else if (Ndim == 3) {
    for (auto k = k_start; k < k_end; k++) {
      for (auto j = j_start; j < j_end; j++) {
        for (auto i = i_start; i < i_end; i++) {
          // result({i, j, k}) = differential_operator(level, i, j, k);
          result({i, j, k}) = ihxx * (solution[level]({i - 1, j, k}) -
                                      2.0 * solution[level]({i, j, k}) +
                                      solution[level]({i + 1, j, k})) +
                              ihyy * (solution[level]({i, j - 1, k}) -
                                      2.0 * solution[level]({i, j, k}) +
                                      solution[level]({i, j + 1, k})) +
                              ihzz * (solution[level]({i, j, k - 1}) -
                                      2.0 * solution[level]({i, j, k}) +
                                      solution[level]({i, j, k + 1}));
        }
      }
    }
  }
}

/**/
template <typename T, int Ndim>
inline void
MultigridBase<T, Ndim>::evaluate_residual(size_t level,
                                          nd::ndArray<T, Ndim> &result) {
  size_t nx = result.size(0);
  size_t ny = result.size(1);
  size_t nz = 0;
  size_t i_start{1}, i_end{nx - 1}, j_start{1}, j_end{ny - 1}, k_start{0},
      k_end{0};
  double hx, hy, hz, ihxx, ihyy, ihzz;

  if (Ndim == 3) {
    nz = result.size(2);
    k_start = 1;
    k_end = nz - 1;
  }
  hx = (lim_end[0] - lim_start[0]) / nx;
  hy = (lim_end[1] - lim_start[1]) / ny;
  hz = (lim_end[2] - lim_start[2]) / nz;
  ihxx = 1.0 / (hx * hx);
  ihyy = 1.0 / (hy * hy);
  ihzz = 1.0 / (hz * hz);

  if (Ndim == 2) {
    for (auto j = j_start; j < j_end; j++) {
      for (auto i = i_start; i < i_end; i++) {
        result({i, j}) =
            source[level]({i, j}) - (ihxx * (solution[level]({i - 1, j}) -
                                             2.0 * solution[level]({i, j}) +
                                             solution[level]({i + 1, j})) +
                                     ihyy * (solution[level]({i, j - 1}) -
                                             2.0 * solution[level]({i, j}) +
                                             solution[level]({i, j + 1})));
      }
    }
  }

  else if (Ndim == 3) {
    for (auto k = k_start; k < k_end; k++) {
      for (auto j = j_start; j < j_end; j++) {
        for (auto i = i_start; i < i_end; i++) {
          result({i, j, k}) = source[level]({i, j, k}) -
                              (ihxx * (solution[level]({i - 1, j, k}) -
                                       2.0 * solution[level]({i, j, k}) +
                                       solution[level]({i + 1, j, k})) +
                               ihyy * (solution[level]({i, j - 1, k}) -
                                       2.0 * solution[level]({i, j, k}) +
                                       solution[level]({i, j + 1, k})) +
                               ihzz * (solution[level]({i, j, k - 1}) -
                                       2.0 * solution[level]({i, j, k}) +
                                       solution[level]({i, j, k + 1})));
        }
      }
    }
  }
}

/**/
// RB-GS smoother  TODO extended it to add JACOBI, SOR, RKL
// Change the interface so it can take as argument an operator Lh <depending of
// the problem both in (2D/3D)>, the unknown vector (e.g defect or potential)
// and the  RHS vector (e.g the source) The name should be change to smoother

template <typename T, int Ndim>
// template <DiscreteOperator Lh>
void MultigridBase<T, Ndim>::smoother(const size_t level, size_t N, T w) {
  size_t nx{solution[level].size(0)}, ny{solution[level].size(1)}, nz{1};
  size_t i_start{1}, i_end{nx - 1}, j_start{1}, j_end{ny - 1}, k_start, k_end;
  double hx, hy, hz, ihxx, ihyy, ihzz;

  if (Ndim == 3) {
    nz = solution[level].size(2);
    k_start = 1;
    k_end = nz - 1;
  }

  hx = (lim_end[0] - lim_start[0]) / nx;
  hy = (lim_end[1] - lim_start[1]) / ny;
  hz = (lim_end[2] - lim_start[2]) / nz;
  ihxx = 1.0 / (hx * hx);
  ihyy = 1.0 / (hy * hy);
  ihzz = 1.0 / (hz * hz);

  // Col-major order
  for (auto round = 0; round < N; round++) {
    for (size_t pColor = 2; pColor > 0;
         pColor--) /* color = 2{RED} color = 1{BLACK}*/
    {
      if (Ndim == 2) {
        for (size_t j = j_start; j < j_end; j++) {
          auto i = pColor - power(-1, pColor) * (j % 2);
          for (; i < i_end; i += 2) {
            relaxation_updater(level, i, j, w);
            // auto sol_tmp =
            //     ((solution[level]({i - 1, j}) + solution[level]({i + 1, j}))
            //     *
            //          ihxx +
            //      (solution[level]({i, j - 1}) + solution[level]({i, j + 1}))
            //      *
            //          ihyy -
            //      source[level]({i, j})) /
            //     (2 * (ihxx + ihyy));
            // solution[level]({i, j}) += w * sol_tmp;
          }
        }
      }

      else if (Ndim == 3) {
        for (size_t k = k_start; k < k_end; k++) {
          for (size_t j = j_start; j < j_end; j++) {
            auto i = pColor - power(-1, pColor) * ((k + j) % 2);
            for (; i < i_end; i += 2) {
              relaxation_updater(level, i, j, k, w);

              //   auto sol_tmp = ((solution[level]({i - 1, j, k}) +
              //                    solution[level]({i + 1, j, k})) *
              //                       ihxx +
              //                   (solution[level]({i, j - 1, k}) +
              //                    solution[level]({i, j + 1, k})) *
              //                       ihyy +
              //                   (solution[level]({i, j, k - 1}) +
              //                    solution[level]({i, j, k + 1})) *
              //                       ihzz -
              //                   source[level]({i, j, k})) /
              //                  (2 * (ihxx + ihyy + ihzz));
              //   solution[level]({i, j, k}) += w * sol_tmp;
            }
          }
        }
      }
    }
    // TODO : Boundaries update
    /** For IMG used zero-fixed boundary conditions at coarser levels.
     *  For FAS used specified physical boundary condition
     */
    if (level != finest_level && bc == ZeroGradient)
      set_zero_gradient_bc<T, Ndim>(solution[level]);
  }
}

/**/
// RB-GS smoother
template <typename T, int Ndim>
// template <DiscreteOperator Lh>
void MultigridBase<T, Ndim>::smoother(const size_t level,
                                      const double tolerance, T w) {
  size_t nx{solution[level].size(0)}, ny{solution[level].size(1)}, nz{1};
  size_t i_start{1}, i_end{nx - 1}, j_start{1}, j_end{ny - 1}, k_start, k_end;

  for (size_t iter = 0; iter < max_iter; iter++) {
    T l2_res = 0, l2_sol = 0, cur_res = 0;
    /*Col-major mode*/

    for (size_t pColor = 2; pColor > 0; pColor--) {
      if (Ndim == 2) {
        for (size_t j = j_start; j < j_end; j++) {
          auto i = pColor - power(-1, pColor) * (j % 2);
          for (; i < i_end; i += 2) {
            cur_res = solution[level]({i, j});
            relaxation_updater(level, i, j, w);

            // current defect
            // TODO include boundaries cells in the residuals norm
            cur_res = solution[level]({i, j}) - cur_res;
            l2_res += power(cur_res, 2);
            l2_sol += power(solution[level]({i, j}), 2);
          }
        }
      }

      else if (Ndim == 3) {
        for (size_t k = k_start; k < k_end; k++) {
          for (size_t j = j_start; j < j_end; j++) {
            auto i = pColor - power(-1, pColor) * ((k + j) % 2);
            for (; i < i_end; i += 2) {
              cur_res = solution[level]({i, j, k});
              relaxation_updater(level, i, j, k, w);
              // current defect
              // TODO include boundaries cells in the residuals norm
              cur_res = solution[level]({i, j, k}) - cur_res;
              l2_res += power(cur_res, 2);
              l2_sol += power(solution[level]({i, j, k}), 2);
            }
          }
        }
      }

      // TODO : Boundaries update
      /** For IMG used zero-fixed boundary conditions at coarser levels.
       *  For FAS used specified physical boundary condition
       */
      if (level != finest_level && bc == ZeroGradient)
        set_zero_gradient_bc<T, Ndim>(solution[level]);
    }

    if ((sqrt(l2_res) / sqrt(l2_sol)) < tolerance)
      return;
  }
}

/**/
template <typename T, int Ndim>
inline void MultigridBase<T, Ndim>::restrict_solution() {
  for (auto level = finest_level;
       level >
       0 /*coarsest_level*/ /*we prefer to use coarsest_level insteed of 0*/;
       level--) {
    solution.coarsen(level, res_ope);
    // updates boundary conditions
    // set_boundary_condition<T,Ndim>(solution[level], ZeroFixed);
  }
}

/**/
template <typename T, int Ndim>
inline void MultigridBase<T, Ndim>::restrict_source() {
  for (auto level = finest_level;
       level >
       0 /*coarsest_level*/ /*we prefer to use coarsest_level insteed of 0*/;
       level--) {
    source.coarsen(level, res_ope);
  }
}

/**/
template <typename T, int Ndim>
inline void MultigridBase<T, Ndim>::set_solution_to_zero(size_t level) {
  solution[level].set_zero();
}

template <typename T, int Ndim>
inline void MultigridBase<T, Ndim>::set_all_solution_to_zero() {
  for (auto level = finest_level; level >= 0; level--)
    set_solution_to_zero(level);
}

template <typename T, int Ndim>
inline void MultigridBase<T, Ndim>::set_source_to_zero(size_t level) {
  source[level].set_zero();
}

template <typename T, int Ndim>
inline void MultigridBase<T, Ndim>::set_all_source_to_zero() {
  for (auto level = finest_level; level >= 0; level--)
    set_source_to_zero(level);
}

/**/
template <typename T, int Ndim>
inline void MultigridBase<T, Ndim>::one_step_down() {
  if (current_level == coarsest_level)
    return;

  // perform Npre smoothing step
  smoother(current_level, npre, w);

  // evaluate the defect/residual in the current stage
  evaluate_residual(current_level, temp[current_level]);

  // restrict the defect to the next level. The restrict defect is store in the
  // next level source nd::ndArray, so it can be used as source for the next
  // level step
  temp.coarsen(current_level, source[current_level - 1], res_ope);

  // set the initial guess for the next level to zero
  set_solution_to_zero(current_level - 1);
  // set boundary condition
  //  set_boundary_condition<T,Ndim>(solution[current_level-1], bc);

  current_level--;
}

/**/

template <typename T, int Ndim>
inline void MultigridBase<T, Ndim>::one_step_up() {
  // prolongate the correction to the finest level
  solution.refine(current_level, temp[current_level + 1], int_ope);

  // apply correction
  solution[current_level + 1] += temp[current_level + 1];

  // perform Npost relaxation
  smoother(current_level + 1, npost, w);
  current_level++;
}

/**/

template <typename T, int Ndim>
inline void MultigridBase<T, Ndim>::mgi_Vcycle() {
  size_t start_level = current_level;
  // Downstroke
  while (current_level != 0 /*coarsest_level*/) {
    one_step_down();
  }
  // solve the coarsest level
  smoother(current_level, res_tol, w);
  // Upstroke
  while (current_level != start_level) {
    one_step_up();
  }
}

/**/
template <typename T, int Ndim>
inline void MultigridBase<T, Ndim>::mgi(const CycleType &cycle_type) {
  if (cycle_type == CycleType::vCycle) {
    mgi_Vcycle();
  } else {
    std::cout << "user-defined multigrid cycle not implemented " << std::endl;
  }

  //   if (cycle_type == CycleType::fCycle) {
  //     mgi_Fcycle();
  //   } else if (cycle_type == CycleType::vCycle) {
  //     mgi_Vcycle();
  //   } else if (cycle_type == CycleType::wCycle) {
  //     mgi_Wcycle();
  //   } else {
  //     std::cout << "user-defined multigrid cycle not implemented " <<
  //     std::endl;
  //   }
}

/**/
template <typename T, int Ndim>
inline void MultigridBase<T, Ndim>::fmg_cycle(
    const CycleType &cycle_type,
    const int mgc_per_level /*Number of multigrid cycle per level*/) {
  for (int level = 0; level < finest_level; level++) {
    TransfertOperator Op = TreeCubic;
    solution.refine(level, Op);
    current_level = level + 1;
    for (int i = 0; i < mgc_per_level; i++)
      mgi(cycle_type);
  }
}

/**/
template <typename T, int Ndim>
inline void MultigridBase<T, Ndim>::fmg(
    const CycleType &cycle_type,
    const int mgc_per_level /*Number of multigrid cycle per level*/) {
  // perform a restriction from the finest level to the coarsest one
  restrict_solution();
  restrict_source();
  // solve the coarsest level
  smoother(coarsest_level, res_tol, w);
  // perform the full multigrid cycle
  fmg_cycle(cycle_type, mgc_per_level);

  /**
   * If in the full multigrid algorithm we apply one V(or W or F)-cycle
   * per level, then at the end of the full multigrid algorithm we need to
   * apply additional V(or W or F)-cycle iterations to obtain a fully-converged
   * solution.
   *
   * See also remark 2.6.2 from Trottenberg et al.2001
   */
  if (fixedNiter)
    iterate_mg_cycle_to_convergence_tol(cycle_type, res_tol);
  else
    iterate_mg_cycle_to_convergence_niter(cycle_type, max_iter);
}
/**/
template <typename T, int Ndim>
inline void MultigridBase<T, Ndim>::iterate_mg_cycle_to_convergence_tol(
    const CycleType &cycle_type, const T eps) {

  nd::ndArray<T, Ndim> _result; /** size of the grid to be specified*/
  int n_iter = 0;

  // computed the residual or the defect ====> should be modified to implement
  // the equation (20) int Tomida & Stone (2023)
  evaluate_residual(finest_level, _result);

  // compute the norm of the residual or the defect

  T def_norm = L2Norm<T, Ndim>(_result);
  while (def_norm > eps) {
    mgi(cycle_type);
    evaluate_residual(finest_level, _result);
    T def_norm_new = L2Norm<T, Ndim>(_result);

    if (def_norm_new / def_norm > 1.0) {
      std::cout
          << "MultigridBase<T,Ndim>::iterateMGCycleToConvergence is diverging"
          << std::endl;
      break;
    }
    if (n_iter > 100) {
      std::cout << "MultigridBase<T,Ndim>::iterateMGCycleToConvergence  take "
                   "too much time : "
                << std::endl;
      std::cout << "Niter done :  " << n_iter << std::endl;
      break;
    }
    n_iter++;
  }
}

/**/
template <typename T, int Ndim>
inline void MultigridBase<T, Ndim>::iterate_mg_cycle_to_convergence_niter(
    const CycleType &cycle_type, const int n_iter) {

  nd::ndArray<T, Ndim> _result; /** size of the grid to be specified*/

  for (int i = 0; i < n_iter; i++) {
    mgi(cycle_type);
    // computed the residual or the defect ====> should be modified to implement
    // the equation (20) int Tomida & Stone (2023)
    evaluate_residual(finest_level, _result);
    T def_norm_new = L2Norm<T, Ndim>(_result);

    // TODO add in the function signature a std::vector<T> to store the history
    // of the residual
  }
}

} // namespace multigrid

#endif