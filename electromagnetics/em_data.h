//
// Created by cepheid on 9/23/24.
//

#ifndef EM_DATA_H
#define EM_DATA_H

#include "core/typelist.h"
#include "aydenstuff/array.h"
#include "core/debug.h"
#include "em_emtpyarray.h"
#include "em_traits.h"
#include "boundaries.h"

using tf::types::Array1D;
using tf::types::Array2D;

template<FieldComponent EXF, FieldComponent EYF, FieldComponent EZF, FieldComponent HXF, FieldComponent HYF, FieldComponent HZF>
struct EMData {
  using value_t = typename EXF::array_t::value_t;
  using dimension_t = typename EXF::array_t::dimension_t;
  using empty_t = EmptyArray<value_t, dimension_t::value>;
  
  using ex_t = typename EXF::array_t;
  using ey_t = typename EYF::array_t;
  using ez_t = typename EZF::array_t;
  using hx_t = typename HXF::array_t;
  using hy_t = typename HYF::array_t;
  using hz_t = typename HZF::array_t;
  
  EMData() = default;
  
  explicit EMData(const size_t nx, const value_t cfl)
  :
    // Ex{}, Jx{}, Cexe{}, Cexh{}, Cjx{},
    // Ey{}, Jy{}, Ceye{}, Ceyh{}, Cjy{},
    Ez{nx}, Jz{nx}, Ceze{nx}, Cezh{nx}, Cjz{nx},
    // Hx{}, Chxe{}, Chxh{},
    Hy{nx - 1}, Chye{nx - 1}, Chyh{nx - 1}
    // Hz{}, Chze{}, Chzh{}
  {
    init_coefficients_1D(cfl);
  }
  
  explicit EMData(const size_t nx, const size_t ny, const value_t cfl)
  :
    // Ex{nx, ny}, Jx{nx, ny}, Cexe{nx, ny}, Cexh{nx, ny}, Cjx{nx, ny},
    // Ey{nx, ny}, Jy{nx, ny}, Ceye{nx, ny}, Ceyh{nx, ny}, Cjy{nx, ny},
    Ez{nx, ny}, Jz{nx, ny}, Ceze{nx, ny}, Cezh{nx, ny}, Cjz{nx, ny},
    Hx{nx, ny - 1}, Chxe{nx, ny - 1}, Chxh{nx, ny - 1},
    Hy{nx - 1, ny}, Chye{nx - 1, ny}, Chyh{nx - 1, ny}
    // Hz{nx, ny}, Chze{nx, ny}, Chzh{nx, ny}
  {
    init_coefficients_2D(cfl);
  }

  void init_coefficients_1D(auto);
  void init_coefficients_2D(auto);

  ex_t Ex;
  ex_t Jx;
  ex_t Cexe;
  ex_t Cexh;
  ex_t Cjx;
  
  ey_t Ey;
  ey_t Jy;
  ey_t Ceye;
  ey_t Ceyh;
  ey_t Cjy;
  
  ez_t Ez;
  ez_t Jz;
  ez_t Ceze;
  ez_t Cezh;
  ez_t Cjz;
  
  hx_t Hx;
  hx_t Chxe;
  hx_t Chxh;
  
  hy_t Hy;
  hy_t Chye;
  hy_t Chyh;
  
  hz_t Hz;
  hz_t Chze;
  hz_t Chzh;
};

template<typename Array>
void init_1d_coeff(Array& arr, auto val) {
  for (size_t i = 0; i < arr.nx; ++i) {
    arr[i] = val;
  }
}

template<typename Array>
requires is_empty_field<Array, EmptyArray1D<typename Array::value_t>>
void init_1d_coeff(Array&, auto) {}

template<typename Array>
void init_2d_coeff(Array& arr, auto val) {
  for (size_t i = 0; i < arr.nx; ++i) {
    for (size_t j = 0; j < arr.nz; ++j) {
      arr(i, j) = val;
    }
  }
}

template<typename Array>
requires is_empty_field<Array, EmptyArray2D<typename Array::value_t>>
void init_2d_coeff(Array&, auto) {}


template <FieldComponent EXF, FieldComponent EYF, FieldComponent EZF, FieldComponent HXF, FieldComponent HYF, FieldComponent HZF>
void EMData<EXF, EYF, EZF, HXF, HYF, HZF>::init_coefficients_1D(auto cfl)
{
  constexpr auto imp0 = 377.0;

  init_1d_coeff(Cexe, 1.0);
  init_1d_coeff(Cexh, cfl * imp0);
  init_1d_coeff(Cjx, 1.0);

  init_1d_coeff(Ceye, 1.0);
  init_1d_coeff(Ceyh, cfl * imp0);
  init_1d_coeff(Cjy, 1.0);

  init_1d_coeff(Ceze, 1.0);
  init_1d_coeff(Cezh, cfl * imp0);
  init_1d_coeff(Cjz, 1.0);

  init_1d_coeff(Chxh, 1.0);
  init_1d_coeff(Chxe, cfl / imp0);

  init_1d_coeff(Chyh, 1.0);
  init_1d_coeff(Chye, cfl / imp0);

  init_1d_coeff(Chzh, 1.0);
  init_1d_coeff(Chze, cfl / imp0);
}

template <FieldComponent EXF, FieldComponent EYF, FieldComponent EZF, FieldComponent HXF, FieldComponent HYF, FieldComponent HZF>
void EMData<EXF, EYF, EZF, HXF, HYF, HZF>::init_coefficients_2D(auto cfl)
{
  constexpr auto imp0 = 377.0;

  init_2d_coeff(Cexe, 1.0);
  init_2d_coeff(Cexh, cfl * imp0);
  init_2d_coeff(Cjx, 1.0);

  init_2d_coeff(Ceye, 1.0);
  init_2d_coeff(Ceyh, cfl * imp0);
  init_2d_coeff(Cjy, 1.0);

  init_2d_coeff(Ceze, 1.0);
  init_2d_coeff(Cezh, cfl * imp0);
  init_2d_coeff(Cjz, 1.0);

  init_2d_coeff(Chxh, 1.0);
  init_2d_coeff(Chxe, cfl / imp0);

  init_2d_coeff(Chyh, 1.0);
  init_2d_coeff(Chye, cfl / imp0);

  init_2d_coeff(Chzh, 1.0);
  init_2d_coeff(Chze, cfl / imp0);
}
#endif //EM_DATA_H
