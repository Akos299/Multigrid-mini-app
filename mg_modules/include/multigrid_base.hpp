#ifndef __MULTIGRID_HPP__
#define __MULTIGRID_HPP__

#include <array>
#include <cstddef>
#include <iostream>

#include "boundary_condition.hpp"
#include "mgrid_stack.hpp"
#include "ndArray.h"
#include "utilities.hpp"

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
    current_level = finest_level;

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
              << " nz_fine = " << nz_fine << "\n";

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
    temp = solution;
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

  /** @brief 
  *
  */
  inline void add_correction(nd::ndArray<T, Ndim> &to_add,int level);

  /** @brief 
  *
  */
  inline void set_solution_at_lev_to_zero(size_t level);

  /** @brief 
  *
  */
  inline void set_source_at_lev_to_zero(size_t level);

  /** @brief 
  *
  */
  inline void set_level_solution(size_t level, T value);
  /** @brief 
  *
  */
  inline void set_level_source(size_t level, T value);
  /** @brief 
  *
  */
  inline void print_level_data(size_t level, std::string data_value);

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
  std::cout << "==== In get_solution === \n ";
  return solution[finest_level];
}

/**/
template <typename T, int Ndim>
template <typename U>
inline void MultigridBase<T, Ndim>::set_initial_guess(U arg) {
  std::cout << "==== In set_initial_guess ==== \n";
  solution[finest_level] = arg;
  is_initial_value_set = true;
}

/**/
template <typename T, int Ndim>
inline nd::ndArray<T, Ndim> &MultigridBase<T, Ndim>::get_source() {
  std::cout << "==== In get_source ==== \n";
  return source[finest_level];
}

/**/
template <typename T, int Ndim>
template <typename U>
inline void MultigridBase<T, Ndim>::set_source(U arg) {
  std::cout << "==== In set_solution ==== \n";
  source[finest_level] = arg;
  is_source_set = true;
}

/**/
template <typename T, int Ndim>
inline void
MultigridBase<T, Ndim>::evaluate_operator(size_t level,
                                          nd::ndArray<T, Ndim> &result) {

  std::cout<< "==== In evaluate_operator ==== \n";
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
  std::cout << "                         === evaluate_residual::start === \n\n";
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
    // std::cout << " buffer infos ( result, source, solution) \n";
    // result.info();
    // source[level].info();
    // solution[level].info();
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
  std::cout << "                         === evaluate_residual::end === \n\n";
}

/**/
// RB-GS smoother  TODO extended it to add JACOBI, SOR, RKL
// Change the interface so it can take as argument an operator Lh <depending of
// the problem both in (2D/3D)>, the unknown vector (e.g defect or potential)
// and the  RHS vector (e.g the source) The name should be change to smoother

template <typename T, int Ndim>
// template <DiscreteOperator Lh>
void MultigridBase<T, Ndim>::smoother(const size_t level, size_t N, T w) {

std::cout << "                        === smoother-with-iter::start_call === \n\n";
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
  std::cout << "                        === smoother-with-iter::end_call === \n\n";
}

/**/
// RB-GS smoother
template <typename T, int Ndim>
// template <DiscreteOperator Lh>
void MultigridBase<T, Ndim>::smoother(const size_t level,
                                      const double tolerance, T w) {
  std::cout<< "  ====@@@@ smoother-with-tol::start_call @@@@==== \n\n";
  size_t nx{solution[level].size(0)}, ny{solution[level].size(1)}, nz{1};
  size_t i_start{1}, i_end{nx - 1}, j_start{1}, j_end{ny - 1}, k_start, k_end;
  if(Ndim == 3)
  {
    nz = solution[level].size(2);
    k_start = 1; k_end = nz - 1;
  }
  int count = 0;
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

 
        // std::cout << " =======start====== \n";
        // std::cout << i_start << " " << j_start << " " << k_start << "\n";
        // std::cout << " =======end======== \n";
        // std::cout << i_end << " " << j_end << " " << k_end << "\n";

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

              // if(count < 10)
              //   std::cout << "smoother::count " << count << "\n";
              // count ++;

            }
          }
        }
      }

      // std::cout << "======================smoother::pcolor============ " << pColor <<'\n';

      // TODO : Boundaries update
      /** For IMG used zero-fixed boundary conditions at coarser levels.
       *  For FAS used specified physical boundary condition
       */
      //  std::cout << "smoother::set_zero_gradient \n";
      if (level != finest_level && bc == ZeroGradient)
        set_zero_gradient_bc<T, Ndim>(solution[level]);
    }
  //  std::cout << "smoother::check_cvg \n";
    if ((sqrt(l2_res) / sqrt(l2_sol)) < tolerance)
      return;
    // std::cout << "smoother::end::check_cvg " << " iter = " << iter << "/" << max_iter << "\n";

  }
   std::cout << "  ====@@@@ smoother-with-tol::end_call   @@@@==== \n\n";
}

/**/
template <typename T, int Ndim>
inline void MultigridBase<T, Ndim>::restrict_solution() {
  std::cout<< " ====@@@@ restrict_solution(all levels)::start @@@@==== \n\n ";
  for (auto level = finest_level;
       level >
       0 /*coarsest_level*/ /*we prefer to use coarsest_level insteed of 0*/;
       level--) {
    solution[level].info();
    solution.coarsen(level, res_ope);
    // updates boundary conditions
    // set_boundary_condition<T,Ndim>(solution[level], ZeroFixed);
  }
  solution[0].info();
  std::cout<< "====@@@@ restrict_solution(all levels)::end @@@@==== \n\n ";
}

/**/
template <typename T, int Ndim>
inline void MultigridBase<T, Ndim>::restrict_source() {
  std::cout<< "  ====@@@@ restrict_source(all levels)::start @@@@==== \n\n ";
  for (auto level = finest_level;
       level >
       0 /*coarsest_level*/ /*we prefer to use coarsest_level insteed of 0*/;
       level--) {
    source.coarsen(level, res_ope);
  }
  std::cout<< " ====@@@@ restrict_source(all levels)::end  @@@@==== \n\n ";
}

/**/
template <typename T, int Ndim>
inline void MultigridBase<T, Ndim>::set_solution_to_zero(size_t level) {
  std::cout << "                         === set_solution_to_zero::level : " << level << "\n\n"; 
  std::cout << "                         === set_solution_to_zero::numel : " << solution[level].numel() << "\n\n";
  std::cout << "                         === set_solution_to_zero::info \n ";
  solution[level].info();
  set_solution_at_lev_to_zero(level);
  std::cout << "                         === set_solution_to_zero::end \n\n"; 
}

template <typename T, int Ndim>
inline void MultigridBase<T, Ndim>::set_all_solution_to_zero() {
  for (auto level = finest_level; level >= 0; level--)
    set_solution_to_zero(level);
}

template <typename T, int Ndim>
inline void MultigridBase<T, Ndim>::set_source_to_zero(size_t level) {
  set_source_at_lev_to_zero(level);
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
  std::cout << "                    ========== one_step_down::smoother ========= \n\n"; 
  smoother(current_level, npre, w);
  std::cout << "                    ========== one_step_down::solution_shape at current level : " << solution[current_level].size(0);
  // evaluate the defect/residual in the current stage
  std::cout << "                    ========== one_step_down::evaluate_residual ==== \n\n";
  evaluate_residual(current_level, temp[current_level]);

  // restrict the defect to the next level. The restrict defect is store in the
  // next level source nd::ndArray, so it can be used as source for the next
  // level step
  std::cout << "                    ========== one_step_down::coarsen =========== \n\n";
  temp.coarsen(current_level, source[current_level - 1], res_ope);


  // set the initial guess for the next level to zero
  std::cout << "                    ========== one_step_down::set_solution_to_zero ==== \n\n";

  /** implement the function to explicitly set solution to zero*/
  set_solution_to_zero(current_level - 1);
  // set boundary condition
  //  set_boundary_condition<T,Ndim>(solution[current_level-1], bc);


  current_level--;
  std::cout << "                   ========== one_step_down::current_level::end ==== \n\n";
}

/**/

template <typename T, int Ndim>
inline void MultigridBase<T, Ndim>::one_step_up() {
  // prolongate the correction to the finest level
  std::cout << "                 ======== one_step_up::refine(interpolation) ==== \n\n";
  solution.refine(current_level, temp[current_level + 1], int_ope);

  temp[current_level+1].info();

  // make a function for the correction step
  // add_correction(solution[current_level + 1], temp[current_level+1]);

  // apply correction
   std::cout << "                ======== one_step_up::correction ====== \n\n";
  // solution[current_level + 1] += temp[current_level + 1];
  add_correction(temp[current_level+1], current_level+1 );

  // perform Npost relaxation
  std::cout << "                 ======= one_step_up::post_smoother ====== \n\n";
  smoother(current_level + 1, npost, w);
  current_level++;
  std::cout << "                 ======= one_step_up::current_level(after one upstroke) ====== :  " <<current_level <<" \n\n";
}

/**/

template <typename T, int Ndim>
inline void MultigridBase<T, Ndim>::mgi_Vcycle() {
  size_t start_level = current_level;
  std::cout << "                 ============ mgi::Vcycle::current_level : " << current_level << "\n\n";
  // Downstroke
  while (current_level != 0 /*coarsest_level*/) {
    one_step_down();
    std::cout << "               =========== mgi_Vcycle::downstroke::cureent_level(after one restriction) :  " << current_level << "\n\n";
  }
  // solve the coarsest level
  
  std::cout << "               ============ mgi_Vcycle::smoother_call (solve the coarsest level)::start =========== \n\n";
  smoother(current_level, res_tol, w);
  std::cout << "               ============ mgi_Vcycle::smoother_call (solve the coarsest level)::end =========== \n\n";
  // Upstroke
  while (current_level != start_level) {
    std::cout << "              ============ mgi_Vcycle::upstroke::cureent_level(start of the upstroke) : " << current_level << "\n\n";
    one_step_up();
    std::cout << "              ============ mgi_Vcycle::upstroke::end ==== \n\n " ;
  }
}

/**/
template <typename T, int Ndim>
inline void MultigridBase<T, Ndim>::mgi(const CycleType &cycle_type) {
  if (cycle_type == CycleType::vCycle) {
    std::cout << "        ==== mgi::vCycle ====     \n\n";  
    mgi_Vcycle();
  } else {
    std::cout << "   ==== user-defined multigrid cycle not implemented ====   \n\n" << std::endl;
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
  std::cout << "  ====@@@@ fmg_cycle::finest_level : " << finest_level << "\n\n";
  for (int level = 0; level < finest_level; level++) {
    TransfertOperator Op = TreeCubic;
    std::cout << " ====@@@@ fmg_cycle::refine at lev (fmg prolongation) : " << level << "\n\n"; 
    solution.refine(level, Op);
    current_level = level + 1;
    std::cout<< " ====@@@@ fmg_cycle::mgc_per_cycle(start mgi iteration) : " << mgc_per_level << "\n\n";
    for (int i = 0; i < mgc_per_level; i++)
      mgi(cycle_type);
  }
}

/**/
template <typename T, int Ndim>
inline void MultigridBase<T, Ndim>::fmg(
  
    const CycleType &cycle_type,
    const int mgc_per_level /*Number of multigrid cycle per level*/) {
  std::cout << "=========================@@@@@ fmg::start ============================@@@@@ \n\n\n\n";
  // perform a restriction from the finest level to the coarsest one
  std::cout << "              =============== fmg::restrict::solution ================= \n\n";
  restrict_solution();
  std::cout << "             =============== fmg::restrict::source    ================ \n\n";
  restrict_source();
  // solve the coarsest level
  std::cout << "             =============== fmg::smoother            ================ \n\n";
  smoother(coarsest_level, res_tol, w);

  // perform the full multigrid cycle
  std::cout << "              ============== fmg::fmg_cycle            ================ \n\n";
  fmg_cycle(cycle_type, mgc_per_level);

  std::cout << "              ============= fmg::begin_iterate         =============== \n\n\n\n";
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
  std::cout << "=========================@@@@@ fmg::end ============================@@@@@ \n\n\n\n";
}
/**/
template <typename T, int Ndim>
inline void MultigridBase<T, Ndim>::iterate_mg_cycle_to_convergence_tol(
    const CycleType &cycle_type, const T eps) {

  std::cout << "==== iterate_mg_cycle_to_convergence_tol::begin ==== \n\n";
  nd::index_t nx{solution[finest_level].size(0)}, ny{solution[finest_level].size(1)}, nz{solution[finest_level].size(2)};
  if (Ndim == 2)
      nz = 1;
  nd::index_t finest_size_3[3] = {nx,ny,nz}, finest_size_2[2]={nx,ny};
  nd::index_t finest_numel   = nx*ny*nz;
  nd::ndArray<T, Ndim> _result;
  if(Ndim == 2)
    _result = nd::ndArray<T, Ndim>(new T[finest_numel], finest_size_2, true) ; /** size of the grid to be specified*/
  else if(Ndim==3)
    _result = nd::ndArray<T, Ndim>(new T[finest_numel], finest_size_3, true);

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
  std::cout << "==== iterate_mg_cycle_to_convergence_tol::end ==== \n";
}

/**/
template <typename T, int Ndim>
inline void MultigridBase<T, Ndim>::iterate_mg_cycle_to_convergence_niter(
    const CycleType &cycle_type, const int n_iter) {
  std::cout << "==== iterate_mg_cycle_to_convergence_niter::begin ==== \n";
   nd::index_t nx{solution[finest_level].size(0)}, ny{solution[finest_level].size(1)}, nz{solution[finest_level].size(2)};
  if (Ndim == 2)
      nz = 1;
  nd::index_t finest_size_3[3] = {nx,ny,nz}, finest_size_2[2]={nx,ny};
  nd::index_t finest_numel   = nx*ny*nz;
  nd::ndArray<T, Ndim> _result;
  if(Ndim == 2)
    _result = nd::ndArray<T, Ndim>(new T[finest_numel], finest_size_2, true) ; /** size of the grid to be specified*/
  else if(Ndim==3)
    _result = nd::ndArray<T, Ndim>(new T[finest_numel], finest_size_3, true); /** size of the grid to be specified*/

  for (int i = 0; i < n_iter; i++) {
    mgi(cycle_type);
    // computed the residual or the defect ====> should be modified to implement
    // the equation (20) int Tomida & Stone (2023)

    std::cout << "iterate_mg_cycle_to_converge_niter:: evaluate_residual \n";
    evaluate_residual(finest_level, _result);

    std::cout << "iterate_mg_cycle_to_converge_niter:: L2 norm \n";
    T def_norm_new = L2Norm<T, Ndim>(_result);

    // TODO add in the function signature a std::vector<T> to store the history
    // of the residual
  }
  std::cout << "==== iterate_mg_cycle_to_convergence_tol::end ==== \n";
}

/**/
template <typename T, int Ndim>
inline void MultigridBase<T, Ndim>::add_correction(nd::ndArray<T, Ndim>& to_add, int level)
{

  std::cout << "                         ==== add_correction::start ==== \n\n";
  size_t nx{solution[level].size(0)}, ny{solution[level].size(1)}, nz{1};
  size_t i_start{1}, i_end{nx - 1}, j_start{1}, j_end{ny - 1}, k_start, k_end;
  if(Ndim == 3)
  {
    nz = solution[level].size(2);
    k_start = 1; k_end = nz - 1;
  }

  if(Ndim == 2)
  {
    for(auto j = j_start; j < j_end; j++)
    {
      for(auto i = i_start; i < i_end; i++)
      {
        solution[level]({i,j}) += to_add({i,j});
      }
    }
  }

  else if (Ndim==3)
  {
    for(auto k = k_start; k < k_end; k++)
    {
      for(auto j = j_start; j < j_end; j++)
      {
        for(auto i = i_start; i < i_end; i++)
        {
          solution[level]({i,j,k}) += to_add({i,j,k});
        }
      }
    }
  }

  std::cout << "                         ==== add_correction::end ==== \n\n";
}

/**/
template<typename T, int Ndim>
inline void MultigridBase<T, Ndim>::set_solution_at_lev_to_zero(size_t level)
{
  size_t nx = solution[level].size(0), ny = solution[level].size(1), nz = 1; 
  if(Ndim == 2)
  {

    for(nd::index_t j = 0; j < ny; j++)
    {
      for(nd::index_t i = 0; i < nx; i++)
      {
         solution[level]({i,j}) = 0;
      }
    }
  }

  if(Ndim == 3)
  {
    nz  = solution[level].size(2);
    for(nd::index_t k = 0; k < nz; k++)
    {
      for(nd::index_t j = 0; j < ny; j++)
        {
          for(nd::index_t i = 0; i < nx; i++)
          {
            solution[level]({i,j,k}) = 0;
          }
        }
    }
  }
}

/**/
template<typename T, int Ndim>
inline void MultigridBase<T, Ndim>::set_source_at_lev_to_zero(size_t level)
{
  size_t nx = source[level].size(0), ny = source[level].size(1), nz = 1; 
  if(Ndim == 2)
  {
    nz  = source[level].size(2);
    for(nd::index_t j = 0; j < ny; j++)
    {
      for(nd::index_t i = 0; i < nx; i++)
      {
         source[level]({i,j}) = 0;
      }
    }
  }

  if(Ndim == 3)
  {
    for(nd::index_t k = 0; k < nz; k++)
    {
      for(nd::index_t j = 0; j < ny; j++)
        {
          for(nd::index_t i = 0; i < nx; i++)
          {
            source[level]({i,j,k}) = 0;
          }
        }
    }
  }
}

/**/
template <typename T, int Ndim>
inline void MultigridBase<T, Ndim>::set_level_solution(size_t level, T value) {
  size_t nx = this->solution[finest_level].size(0), ny = this->solution[finest_level].size(1), nz = 1;
  if (Ndim == 3) {

    std::cout << "==================================    multigridBase::set_level_solution  ====================== \n";
    nz = this->solution[finest_level].size(2);

    std::cout << " nx " << nx << " ny " << ny << " nz " << nz << "\n";
    for (nd::index_t k = 0; k < nz; k++) {
      for (nd::index_t j = 0; j < ny; j++) {
        for (nd::index_t i = 0; i < nx; i++) {
          solution[finest_level]({i, j, k}) = value;
        }
      }
    }
  }

  if (Ndim == 2) {
    for (nd::index_t j = 0; j < ny; j++) {
      for (nd::index_t i = 0; i < nx; i++) {
        solution[finest_level]({i, j}) = value;
      }
    }
  }
}

/**/
template <typename T, int Ndim>
inline void MultigridBase<T, Ndim>::set_level_source(size_t level, T value) {
  size_t nx = source[finest_level].size(0), ny = source[1].size(1), nz = 1;
  if (Ndim == 3) {
    nz = source[finest_level].size(2);
    for (nd::index_t k = 0; k < nz; k++) {
      for (nd::index_t j = 0; j < ny; j++) {
        for (nd::index_t i = 0; i < nx; i++) {
          source[finest_level]({i, j, k}) = value;
        }
      }
    }
  }

  if (Ndim == 2) {
    for (nd::index_t j = 0; j < ny; j++) {
      for (nd::index_t i = 0; i < nx; i++) {
        source[finest_level]({i, j}) = value;
      }
    }
  }
}


/**/
template <typename T, int Ndim>
inline void
MultigridBase<T, Ndim>::print_level_data(size_t level,
                                             std::string data_value) {
  size_t n_elts = 0;
  if (data_value == "solution") {
    for (nd::index_t k = 0; k < solution[level].size(2); k++) {
      for (nd::index_t j = 0; j < solution[level].size(1); j++) {
        for (nd::index_t i = 0; i < solution[level].size(0); i++) {
          std::cout << solution[level]({i, j, k}) << " ";
        }
      }
    }
    std::cout << "\n";
  }
  if (data_value == "source") {
    for (nd::index_t k = 0; k < source[level].size(2); k++) {
      for (nd::index_t j = 0; j < source[level].size(1); j++) {
        for (nd::index_t i = 0; i < source[level].size(0); i++) {
          std::cout << source[level]({i, j, k}) << " ";
        }
      }
    }
    std::cout << "\n";
  }
}

} // namespace multigrid

#endif