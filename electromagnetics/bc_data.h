//
// Created by cepheid on 10/23/24.
//

#ifndef BC_DATA_H
#define BC_DATA_H

#include "em_traits.h"
#include "aydenstuff/array.h"

using tf::types::Array1D;
using tf::types::Array2D;
using tf::types::Array3D;

template<typename array_t>
struct PMLData {
  using value_t = typename array_t::value_t;
  using dimension_t = typename array_t::dimension_t;
  
  explicit PMLData(const size_t nx) requires dimension_t::value == 1
  : psi{nx}, b{nx}, c{nx}
  {}

  explicit PMLData(const size_t nx, const size_t ny) requires dimension_t::value == 2
  : psi{nx, ny}, b{nx, ny}, c{nx, ny}
  {}

  explicit PMLData(const size_t nx, const size_t ny, const size_t nz) requires dimension_t::value == 3
  : psi{nx - 1, ny, nz}, b{nx - 1, ny, nz}, c{nx - 1, ny, nz}
  {}
  
  array_t psi;
  array_t b;
  array_t c;
};
  
template<typename T, FieldComponent EXF, FieldComponent EYF, FieldComponent EZF, FieldComponent HXF, FieldComponent HYF, FieldComponent HZF>
struct BCData {
  using array_t = T;
  using value_t = typename array_t::value_t;
  using dimension_t = typename array_t::dimension_t;
  using empty_t = EmptyArray<value_t, dimension_t::value>;
  
  using ex_t = PMLData<typename EXF::array_t>;
  using ey_t = PMLData<typename EYF::array_t>;
  using ez_t = PMLData<typename EZF::array_t>;
  using hx_t = PMLData<typename HXF::array_t>;
  using hy_t = PMLData<typename HYF::array_t>;
  using hz_t = PMLData<typename HZF::array_t>;

  explicit BCData(const size_t nx) requires (dimension_t::value == 1)
  // : psiEz{nx}, bEz{nx}, cEz{nx},
  //   psiHy{nx - 1}, bHy{nx - 1}, cHy{nx - 1}
  {
    // init_coefficients_1D(nx, cfl);
  }

  explicit BCData(const size_t nx, const size_t ny) requires (dimension_t::value == 2 and !is_empty_field<hx_t, empty_t>)
  // : psiEz{nx, ny}, bEz{nx, ny}, cEz{nx, ny},
  //   psiHx{nx, ny - 1}, bHx{nx, ny - 1}, cHx{nx, ny - 1},
  //   psiHy{nx - 1, ny}, bHy{nx - 1, ny}, cHy{nx - 1, ny}
  {
    // init_coefficients_2D(nx, cfl);
  }

  explicit BCData(const size_t nx, const size_t ny, const size_t nz) requires (dimension_t::value == 3)
  // : psiEx{nx - 1, ny, nz}, bEx{nx - 1, ny, nz}, cEx{nx - 1, ny, nz},
  //   psiEy{nx, ny - 1, nz}, bEy{nx, ny - 1, nz}, cEy{nx, ny - 1, nz},
  //   psiEz{nx, ny, nz - 1}, bEz{nx, ny, nz - 1}, cEz{nx, ny, nz - 1},
  //   psiHx{nx, ny - 1, nz - 1}, bHx{nx, ny - 1, nz - 1}, cHx{nx, ny - 1, nz - 1},
  //   psiHy{nx - 1, ny, nz - 1}, bHy{nx - 1, ny, nz - 1}, cHy{nx - 1, ny, nz - 1},
  //   psiHz{nx - 1, ny - 1, nz}, bHz{nx - 1, ny - 1, nz}, cHz{nx - 1, ny - 1, nz}
  {
    // init_coefficients_2D(nx, cfl);
  }

  ex_t Ex_pml_Y0;
  ex_t Ex_pml_Y1;
  ex_t Ex_pml_Z0;
  ex_t Ex_pml_Z1;

  ey_t Ey_pml_X0;
  ey_t Ey_pml_X1;
  ey_t Ey_pml_Z0;
  ey_t Ey_pml_Z1;

  ex_t Ez_pml_X0;
  ex_t Ez_pml_X1;
  ex_t Ez_pml_Y0;
  ex_t Ez_pml_Y1;
  
  // ex_t psiEx;
  // ex_t bEx;
  // ex_t cEx;
  //
  // ey_t psiEy;
  // ey_t bEy;
  // ey_t cEy;
  //
  // ez_t psiEz;
  // ez_t bEz;
  // ez_t cEz;
  //
  // hx_t psiHx;
  // hx_t bHx;
  // hx_t cHx;
  //
  // hy_t psiHy;
  // hy_t bHy;
  // hy_t cHy;
  //
  // hz_t psiHz;
  // hz_t bHz;
  // hz_t cHz;
};

#endif //BC_DATA_H
