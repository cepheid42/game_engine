//
// Created by cepheid on 9/23/24.
//

#ifndef EM_DATA_H
#define EM_DATA_H

#include "core/typelist.h"
#include "aydenstuff/array.h"
#include "em_emtpyarray.h"
#include "em_traits.h"

using tf::types::Array1D;
using tf::types::Array2D;

enum class BoundaryType { None, Periodic, PML };

constexpr BoundaryType BCX0 = BoundaryType::PML;
constexpr BoundaryType BCY0 = BoundaryType::PML;
constexpr BoundaryType BCZ0 = BoundaryType::PML;
constexpr BoundaryType BCX1 = BoundaryType::PML;
constexpr BoundaryType BCY1 = BoundaryType::PML;
constexpr BoundaryType BCZ1 = BoundaryType::PML;

template<typename T>
std::vector<T> linspace(T start, T stop, size_t n_points, const bool endpoint=true) {
  std::vector<T> result(n_points);
  if (endpoint) {
    n_points -= 1;
    result[result.size() - 1] = stop;
  }
  auto delta = (stop - start) / static_cast<T>(n_points);
  T val = start;
  for (size_t i = 0; i < n_points; ++i) {
    result[i] = val;
    val += delta;
  }
  return result;
}

template<FieldComponent EXF, FieldComponent EYF, FieldComponent EZF, FieldComponent HXF, FieldComponent HYF, FieldComponent HZF>
struct BCData {
  using value_t = typename EXF::value_t;
  using dimension_t = typename EXF::dimension_t;

  using ex_t = typename EXF::array_t;
  using ey_t = typename EYF::array_t;
  using ez_t = typename EZF::array_t;
  using hx_t = typename HXF::array_t;
  using hy_t = typename HYF::array_t;
  using hz_t = typename HZF::array_t;

  explicit BCData(const size_t nx)
  : psiEz{nx}, bEz{nx}, cEz{nx},
    psiHy{nx - 1}, bHy{nx - 1}, cHy{nx - 1}
  {}

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

template<FieldComponent EXF, FieldComponent EYF, FieldComponent EZF, FieldComponent HXF, FieldComponent HYF, FieldComponent HZF>
struct EMData {
  using value_t = typename EXF::array_t::value_t;
  using dimension_t = typename EXF::array_t::dimension_t;

  using boundary_t = BCData<EXF, EYF, EZF, HXF, HYF, HZF>;
  
  using empty_t = EmptyArray<value_t, dimension_t::value>;
  
  using ex_t = typename EXF::array_t;
  using ey_t = typename EYF::array_t;
  using ez_t = typename EZF::array_t;
  using hx_t = typename HXF::array_t;
  using hy_t = typename HYF::array_t;
  using hz_t = typename HZF::array_t;
  
  EMData() = default;
  
  explicit EMData(const size_t nx, const value_t cfl)
  : Ez{nx}, Jz{nx}, Ceze{nx}, Cezh{nx}, Cjz{nx},
    Hy{nx - 1}, Chye{nx - 1}, Chyh{nx - 1},
    lo{nx}, hi{nx}
  {
    init_coefficients_1D(nx, cfl);
  }
  
  // explicit EMData(const size_t nx, const size_t ny, const value_t cfl)
  // :
  //   // Ex{nx, ny}, Jx{nx, ny}, Cexe{nx, ny}, Cexh{nx, ny}, Cjx{nx, ny},
  //   // Ey{nx, ny}, Jy{nx, ny}, Ceye{nx, ny}, Ceyh{nx, ny}, Cjy{nx, ny},
  //   Ez{nx, ny}, Jz{nx, ny}, Ceze{nx, ny}, Cezh{nx, ny}, Cjz{nx, ny},
  //   Hx{nx, ny - 1}, Chxe{nx, ny - 1}, Chxh{nx, ny - 1},
  //   Hy{nx - 1, ny}, Chye{nx - 1, ny}, Chyh{nx - 1, ny}
  //   // Hz{nx, ny}, Chze{nx, ny}, Chzh{nx, ny}
  // {
  //   init_coefficients_2D(cfl);
  // }

  void init_coefficients_1D(auto, auto);
  // void init_coefficients_2D(auto);

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

  boundary_t lo;
  boundary_t hi;  
};

template<typename Array>
requires is_empty_field<Array, EmptyArray1D<typename Array::value_t>>
void init_1d_coeff(Array&, auto) {}

template<typename Array>
void init_1d_coeff(Array& arr, auto val) {
  for (size_t i = 0; i < arr.nx; ++i) {
    arr[i] = val;
  }
}

template<typename Array, bool Upper, BoundaryType Boundary>
requires is_empty_field<Array, EmptyArray1D<typename Array::value_t>> or (Boundary != BoundaryType::PML)
void init_1d_pml(Array&, Array&, auto) {}

template<typename Array, bool Upper, BoundaryType Boundary>
void init_1d_pml(Array& b, Array& c, const auto dx) {
  constexpr size_t dPML = 10u;

  if constexpr (Upper) {
    std::ranges::reverse();
  }
}

// template<typename Array>
// void init_2d_coeff(Array& arr, auto val) {
//   for (size_t i = 0; i < arr.nx; ++i) {
//     for (size_t j = 0; j < arr.nz; ++j) {
//       arr(i, j) = val;
//     }
//   }
// }
//
// template<typename Array>
// requires is_empty_field<Array, EmptyArray2D<typename Array::value_t>>
// void init_2d_coeff(Array&, auto) {}




template <FieldComponent EXF, FieldComponent EYF, FieldComponent EZF, FieldComponent HXF, FieldComponent HYF, FieldComponent HZF>
void EMData<EXF, EYF, EZF, HXF, HYF, HZF>::init_coefficients_1D(auto nx, auto cfl)
{
  constexpr auto imp0 = 377.0;
  constexpr auto dx = 1.0 / static_cast<decltype(cfl)>(nx - 1);

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

  init_1d_pml(lo.bEx, lo.cEx, dx);
  init_1d_pml(lo.bEy, lo.cEy, dx);
  init_1d_pml<EZF::array_t, false, BCX0>(lo.bEz, lo.cEz, dx);
  init_1d_pml(lo.bHx, lo.cHx, dx);
  init_1d_pml<HYF::array_t, false, BCX0>(lo.bHy, lo.cHy, dx);
  init_1d_pml(lo.bHz, lo.cHz, dx);

  init_1d_pml(hi.bEx, hi.cEx, dx);
  init_1d_pml(hi.bEy, hi.cEy, dx);
  init_1d_pml<EZF::array_t, true, BCX0>(hi.bEz, hi.cEz, dx);
  init_1d_pml(hi.bHx, hi.cHx, dx);
  init_1d_pml<HYF::array_t, true, BCX0>(hi.bHy, hi.cHy, dx);
  init_1d_pml(hi.bHz, hi.cHz, dx);
}

// template <FieldComponent EXF, FieldComponent EYF, FieldComponent EZF, FieldComponent HXF, FieldComponent HYF, FieldComponent HZF>
// void EMData<EXF, EYF, EZF, HXF, HYF, HZF>::init_coefficients_2D(auto cfl)
// {
//   constexpr auto imp0 = 377.0;
//
//   init_2d_coeff(Cexe, 1.0);
//   init_2d_coeff(Cexh, cfl * imp0);
//   init_2d_coeff(Cjx, 1.0);
//
//   init_2d_coeff(Ceye, 1.0);
//   init_2d_coeff(Ceyh, cfl * imp0);
//   init_2d_coeff(Cjy, 1.0);
//
//   init_2d_coeff(Ceze, 1.0);
//   init_2d_coeff(Cezh, cfl * imp0);
//   init_2d_coeff(Cjz, 1.0);
//
//   init_2d_coeff(Chxh, 1.0);
//   init_2d_coeff(Chxe, cfl / imp0);
//
//   init_2d_coeff(Chyh, 1.0);
//   init_2d_coeff(Chye, cfl / imp0);
//
//   init_2d_coeff(Chzh, 1.0);
//   init_2d_coeff(Chze, cfl / imp0);
// }

#endif //EM_DATA_H
