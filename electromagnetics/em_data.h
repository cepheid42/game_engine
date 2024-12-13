//
// Created by cepheid on 9/23/24.
//

#ifndef EM_DATA_H
#define EM_DATA_H

#include "aydenstuff/array.h"
#include "em_traits.h"

// todo: Need to add namespaces here

using tf::types::Array1D;
using tf::types::Array2D;
using tf::types::Array3D;

template<FieldComponent ex_t, FieldComponent ey_t, FieldComponent ez_t, FieldComponent hx_t, FieldComponent hy_t, FieldComponent hz_t>
struct EMData {
  using value_t = typename ex_t::value_t;
  using dimension_t = typename ex_t::dimension_t;
  using empty_t = EmptyArray<value_t, dimension_t::value>;

  EMData() = default;

  explicit EMData(const size_t nx, const value_t dt)
  requires (dimension_t::value == 1)
  : Ez{nx}, Jz{nx}, Ceze{nx}, Cezhy{nx}, Cjz{nx},
    Hy{nx - 1}, Chyez{nx - 1}, Chyh{nx - 1}
  {
    init_coefficients(dt);
  }

  explicit EMData(const size_t nx, const size_t ny, const value_t dt)
  requires (dimension_t::value == 2 and !is_empty_field<hx_t, empty_t>) // todo: make this smarter so other directions can be used.
  : Ez{nx, ny}, Jz{nx, ny}, Ceze{nx, ny}, Cezhx{nx, ny}, Cezhy{nx, ny}, Cjz{nx, ny},
    Hx{nx, ny - 1}, Chxez{nx, ny - 1}, Chxh{nx, ny - 1},
    Hy{nx - 1, ny}, Chyez{nx - 1, ny}, Chyh{nx - 1, ny}
  {
    // TMz constructor
    init_coefficients(dt);
  }

  explicit EMData(const size_t nx, const size_t ny, const value_t dt)
  requires (dimension_t::value == 2 and !is_empty_field<ex_t, empty_t>)// todo: make this smarter so other directions can be used.
  : Ex{nx - 1, ny}, Jx{nx - 1, ny}, Cexe{nx - 1, ny}, Cexhz{nx - 1, ny}, Cjx{nx - 1, ny},
    Ey{nx, ny - 1}, Jy{nx, ny - 1}, Ceye{nx, ny - 1}, Ceyhz{nx, ny - 1}, Cjy{nx, ny - 1},
    Hz{nx - 1, ny - 1}, Chzex{nx - 1, ny - 1}, Chzey{nx - 1, ny - 1}, Chzh{nx - 1, ny - 1}
  {
    // TEz constructor
    init_coefficients(dt);
  }

  explicit EMData(const size_t nx, const size_t ny, const size_t nz, const value_t dt)
  requires (dimension_t::value == 3)
  : Ex{nx - 1, ny, nz}, Jx{nx - 1, ny, nz}, Cexe{nx - 1, ny, nz}, Cexhy{nx - 1, ny, nz}, Cexhz{nx - 1, ny, nz}, Cjx{nx - 1, ny, nz},
    Ey{nx, ny - 1, nz}, Jy{nx, ny - 1, nz}, Ceye{nx, ny - 1, nz}, Ceyhx{nx, ny - 1, nz}, Ceyhz{nx, ny - 1, nz}, Cjy{nx, ny - 1, nz},
    Ez{nx, ny, nz - 1}, Jz{nx, ny, nz - 1}, Ceze{nx, ny, nz - 1}, Cezhx{nx, ny, nz - 1}, Cezhy{nx, ny, nz - 1}, Cjz{nx, ny, nz - 1},
    Hx{nx, ny - 1, nz - 1}, Chxey{nx, ny - 1, nz - 1}, Chxez{nx, ny - 1, nz - 1}, Chxh{nx, ny - 1, nz - 1},
    Hy{nx - 1, ny, nz - 1}, Chyex{nx - 1, ny, nz - 1}, Chyez{nx - 1, ny, nz - 1}, Chyh{nx - 1, ny, nz - 1},
    Hz{nx - 1, ny - 1, nz}, Chzex{nx - 1, ny - 1, nz}, Chzey{nx - 1, ny - 1, nz}, Chzh{nx - 1, ny - 1, nz}
  {
    init_coefficients(dt);
  }

  void init_coefficients(value_t);
  // void init_coefficients_2D(value_t);

  ex_t Ex;
  ex_t Jx;
  ex_t Cexe;
  ex_t Cexhy;
  ex_t Cexhz;
  ex_t Cjx;

  ey_t Ey;
  ey_t Jy;
  ey_t Ceye;
  ey_t Ceyhx;
  ey_t Ceyhz;
  ey_t Cjy;

  ez_t Ez;
  ez_t Jz;
  ez_t Ceze;
  ez_t Cezhx;
  ez_t Cezhy;
  ez_t Cjz;

  hx_t Hx;
  hx_t Chxey;
  hx_t Chxez;
  hx_t Chxh;

  hy_t Hy;
  hy_t Chyex;
  hy_t Chyez;
  hy_t Chyh;

  hz_t Hz;
  hz_t Chzex;
  hz_t Chzey;
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
void EMData<EXF, EYF, EZF, HXF, HYF, HZF>::init_coefficients(const value_t dt)
{
  constexpr auto eps0 = 8.854187812813e-12;
  constexpr auto mu0 = 1.2566370621219e-6;
  constexpr auto sigma = 0.0;

  constexpr auto dx = 1.0 / 99.0;

  // half dt for H field, since it's split into two steps
  const auto hc = 0.5 * dt / (mu0 * dx);

  const auto e_num = dt / (eps0 * dx);
  const auto alpha = (sigma * dt) / (2.0 * eps0);
  const auto ec = (1.0 - alpha) / (1.0 + alpha);
  const auto eh = e_num / (1.0 + alpha);

  init_coeff(Cexe, ec);
  init_coeff(Cexhy, eh);
  init_coeff(Cexhz, eh);
  init_coeff(Cjx, dt / eps0);

  init_coeff(Ceye, ec);
  init_coeff(Ceyhx, eh);
  init_coeff(Ceyhz, eh);
  init_coeff(Cjy, dt / eps0);

  init_coeff(Ceze, ec);
  init_coeff(Cezhx, eh);
  init_coeff(Cezhy, eh);
  init_coeff(Cjz, dt / eps0);

  init_coeff(Chxh, 1.0);
  init_coeff(Chxey, hc);
  init_coeff(Chxez, hc);

  init_coeff(Chyh, 1.0);
  init_coeff(Chyex, hc);
  init_coeff(Chyez, hc);

  init_coeff(Chzh, 1.0);
  init_coeff(Chzex, hc);
  init_coeff(Chzey, hc);
}

#endif //EM_DATA_H
