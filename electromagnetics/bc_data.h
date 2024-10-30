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

template<typename T, FieldComponent EXF, FieldComponent EYF, FieldComponent EZF, FieldComponent HXF, FieldComponent HYF, FieldComponent HZF>
struct BCData {
  using array_t = T;
  using value_t = typename array_t::value_t;
  using dimension_t = typename array_t::dimension_t;
  using empty_t = EmptyArray<value_t, dimension_t::value>;
  
  using ex_t = typename EXF::array_t;
  using ey_t = typename EYF::array_t;
  using ez_t = typename EZF::array_t;
  using hx_t = typename HXF::array_t;
  using hy_t = typename HYF::array_t;
  using hz_t = typename HZF::array_t;

  explicit BCData(const size_t nx) requires (dimension_t::value == 1)
  : psiEz{nx}, bEz{nx}, cEz{nx},
    psiHy{nx - 1}, bHy{nx - 1}, cHy{nx - 1}
  {
    // init_coefficients_1D(nx, cfl);
  }

  explicit BCData(const size_t nx, const size_t ny) requires (dimension_t::value == 2 and !is_empty_field<hx_t, empty_t>)
  : psiEz{nx, ny}, bEz{nx, ny}, cEz{nx, ny},
    psiHx{nx, ny - 1}, bHx{nx, ny - 1}, cHx{nx, ny - 1},
    psiHy{nx - 1, ny}, bHy{nx - 1, ny}, cHy{nx - 1, ny}
  {
    // init_coefficients_2D(nx, cfl);
  }

  ex_t psiEx;
  ex_t bEx;
  ex_t cEx;

  ey_t psiEy;
  ey_t bEy;
  ey_t cEy;

  ez_t psiEz;
  ez_t bEz;
  ez_t cEz;

  hx_t psiHx;
  hx_t bHx;
  hx_t cHx;

  hy_t psiHy;
  hy_t bHy;
  hy_t cHy;

  hz_t psiHz;
  hz_t bHz;
  hz_t cHz;
};

#endif //BC_DATA_H
