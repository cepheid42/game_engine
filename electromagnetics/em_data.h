//
// Created by cepheid on 9/23/24.
//

#ifndef EM_DATA_H
#define EM_DATA_H

#include "core/debug.h"
// #include "core/typelist.h"
#include "aydenstuff/array.h"
#include "em_emtpyarray.h"
#include "em_traits.h"
// #include "bc_data.h"

using tf::types::Array1D;
using tf::types::Array2D;
using tf::types::Array3D;

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
  requires (dimension_t::value == 1)
  : Ez{nx}, Jz{nx}, Ceze{nx}, Cezh{nx}, Cjz{nx},
    Hy{nx - 1}, Chye{nx - 1}, Chyh{nx - 1}
  {
    init_coefficients(cfl);
  }
  
  explicit EMData(const size_t nx, const size_t ny, const value_t cfl)
  requires (dimension_t::value == 2 and !is_empty_field<hx_t, empty_t>)
  : Ez{nx, ny}, Jz{nx, ny}, Ceze{nx, ny}, Cezh{nx, ny}, Cjz{nx, ny},
    Hx{nx, ny - 1}, Chxe{nx, ny - 1}, Chxh{nx, ny - 1},
    Hy{nx - 1, ny}, Chye{nx - 1, ny}, Chyh{nx - 1, ny}
  {
    // TMz constructor
    init_coefficients(cfl);
  }

  // explicit EMData(const size_t nx, const size_t ny, const value_t cfl)
  // requires (dimension_t::value == 2 and !is_empty_field<ex_t, empty_t>)
  // : Ez{nx, ny}, Jz{nx, ny}, Ceze{nx, ny}, Cezh{nx, ny}, Cjz{nx, ny}, // todo: wrong fields here
  //   Hx{nx, ny - 1}, Chxe{nx, ny - 1}, Chxh{nx, ny - 1},
  //   Hy{nx - 1, ny}, Chye{nx - 1, ny}, Chyh{nx - 1, ny}
  // {
  //   // TEz constructor
  //   init_coefficients(cfl);
  // }

  explicit EMData(const size_t nx, const size_t ny, const size_t nz, const value_t cfl)
  requires (dimension_t::value == 3)
  : Ex{nx - 1, ny, nz}, Jx{nx - 1, ny, nz}, Cexe{nx - 1, ny, nz}, Cexh{nx - 1, ny, nz}, Cjx{nx - 1, ny, nz},
    Ey{nx, ny - 1, nz}, Jy{nx, ny - 1, nz}, Ceye{nx, ny - 1, nz}, Ceyh{nx, ny - 1, nz}, Cjy{nx, ny - 1, nz},
    Ez{nx, ny, nz - 1}, Jz{nx, ny, nz - 1}, Ceze{nx, ny, nz - 1}, Cezh{nx, ny, nz - 1}, Cjz{nx, ny, nz - 1},
    Hx{nx, ny - 1, nz - 1}, Chxe{nx, ny - 1, nz - 1}, Chxh{nx, ny - 1, nz - 1},
    Hy{nx - 1, ny, nz - 1}, Chye{nx - 1, ny, nz - 1}, Chyh{nx - 1, ny, nz - 1},
    Hz{nx - 1, ny - 1, nz}, Chze{nx - 1, ny - 1, nz}, Chzh{nx - 1, ny - 1, nz}
  {
    init_coefficients(cfl);
  }

  void init_coefficients(value_t);
  // void init_coefficients_2D(value_t);

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
requires is_empty_field<Array, EmptyArray<typename Array::value_t, Array::dimension_t::value>>
void init_coeff(Array&, auto) {}

template<typename Array>
void init_coeff(Array& arr, auto val) {
  for (auto& el: arr) {
    el = val;
  }
}

template <FieldComponent EXF, FieldComponent EYF, FieldComponent EZF, FieldComponent HXF, FieldComponent HYF, FieldComponent HZF>
void EMData<EXF, EYF, EZF, HXF, HYF, HZF>::init_coefficients(const value_t cfl)
{
  constexpr auto imp0 = 377.0;

  init_coeff(Cexe, 1.0);
  init_coeff(Cexh, cfl * imp0);
  init_coeff(Cjx, 1.0);

  init_coeff(Ceye, 1.0);
  init_coeff(Ceyh, cfl * imp0);
  init_coeff(Cjy, 1.0);

  init_coeff(Ceze, 1.0);
  init_coeff(Cezh, cfl * imp0);
  init_coeff(Cjz, 1.0);

  init_coeff(Chxh, 1.0);
  init_coeff(Chxe, cfl / imp0);

  init_coeff(Chyh, 1.0);
  init_coeff(Chye, cfl / imp0);

  init_coeff(Chzh, 1.0);
  init_coeff(Chze, cfl / imp0);
}

#endif //EM_DATA_H
