#ifndef __BCONDITIONS_HPP__
#define __BCONDITIONS_HPP__

#include "../include/utilities.hpp"
#include "../ndarray/ndArray.hpp"

namespace multigrid {

template <typename T, nd::dimen_t Ndim>
inline void set_periodic_bc(nd::ndArray<T, Ndim> &buffer) {
  auto nx{buffer.size(0)}, ny{buffer.size(1)}, nz{0};

  if (Ndim == 2) {
    // set at j = 0, ny-1
    for (int i = 0; i < nx; i++)
      buffer({i, 0}) = buffer({i, ny - 1});
    // set tat i = 0, nx - 1
    for (int j = 0; j < ny; j++)
      buffer({0, j}) = buffer({nx - 1, j});
  } else if (Ndim == 3) {
    // set at j = 0,ny-1
    nz = buffer.size(2);
    for (int k = 0; k < nz; k++)
      for (int i = 0; i < nx; i++)
        buffer({i, 0, k}) = buffer({i, ny - 1, k});
    // set at i = 0, nx-1
    for (int k = 0; k < nz; k++)
      for (int j = 0; j < ny; j++)
        buffer({0, j, k}) = buffer({nx - 1, j, k});
    // set at k = 0, nz-1
    for (int j = 0; j < ny; j++)
      for (int i = 0; i < nx; i++)
        buffer({i, j, 0}) = buffer({i, j, nz - 1});
  }
}

template <typename T, nd::dimen_t Ndim>
inline void set_zero_gradient_bc(nd::ndArray<T, Ndim> &buffer) {
  auto nx{buffer.size(0)}, ny{buffer.size(1)}, nz{0};
  if (Ndim == 2) {
    // set at j = 0, ny-1
    for (int i = 0; i < nx; i++) {
      buffer({i, 0}) = buffer({i, 1});
      buffer({i, ny - 1}) = buffer({i, ny - 2});
    }
    // set tat i = 0, nx - 1
    for (int j = 0; j < ny; j++) {
      buffer({0, j}) = buffer({1, j});
      buffer({nx - 1, j}) = buffer({nx - 2, j});
    }
  } else if (Ndim == 3) {

    nz = buffer.size(2);
    // set at j = 0,ny-1
    for (int k = 0; k < nz; k++) {
      for (int i = 0; i < nx; i++) {
        buffer({i, 0, k}) = buffer({i, 1, k});
        buffer({i, ny - 1, k}) = buffer({i, ny - 2, k});
      }
    }
    // set at i = 0, nx-1
    for (int k = 0; k < nz; k++) {
      for (int j = 0; j < ny; j++) {
        buffer({0, j, k}) = buffer({1, j, k});
        buffer({nx - 1, j, k}) = buffer({nx - 2, j, k});
      }
    }
    // set at k = 0, nz-1
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++) {
        buffer({i, j, 0}) = buffer({i, j, 1});
        buffer({i, j, nz - 1}) = buffer({i, j, nz - 2});
      }
    }
  }
}

template <typename T, nd::dimen_t Ndim>
inline void set_zero_fixed_bc(nd::ndArray<T, Ndim> &buffer) {
  auto nx{buffer.size(0)}, ny{buffer.size(1)}, nz{0};
  if (Ndim == 2) {
    // set at j = 0, ny-1
    for (int i = 0; i < nx; i++) {
      buffer({i, 0}) = 0;
      buffer({i, ny - 1}) = 0;
    }
    // set tat i = 0, nx - 1
    for (int j = 0; j < ny; j++) {
      buffer({0, j}) = 0;
      buffer({nx - 1, j}) = 0;
    }
  } else if (Ndim == 3) {
    nz = buffer.size(2);
    // set at j = 0,ny-1
    for (int k = 0; k < nz; k++) {
      for (int i = 0; i < nx; i++) {
        buffer({i, 0, k}) = 0;
        buffer({i, ny - 1, k}) = 0;
      }
    }
    // set at i = 0, nx-1
    for (int k = 0; k < nz; k++) {
      for (int j = 0; j < ny; j++) {
        buffer({0, j, k}) = 0;
        buffer({nx - 1, j, k}) = 0;
      }
    }
    // set at k = 0, nz-1
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++) {
        buffer({i, j, 0}) = 0;
        buffer({i, j, nz - 1}) = 0;
      }
    }
  }
}

template <typename T, nd::dimen_t Ndim>
inline void set_boundary_condition(nd::ndArray<T, Ndim>& buffer,BcType& bc_type){
    if (bc_type == BcType::Periodic)
        set_periodic_bc(buffer);
    else if(bc_type == BcType::ZeroFixed)
        set_zero_fixed_bc(buffer);
    else if(bc_type == BcType::ZeroGradient)
        set_zero_gradient_bc(buffer);
    else{/** TO BE DEFINED*/}
}
} // namespace multigrid

#endif