//
// Created by cepheid on 6/28/24.
//

#ifndef EM_SOLVER_H
#define EM_SOLVER_H

#include <cassert>
#include <array>
// #include <iostream>
// #include <string_view>

#include "../aydenstuff/array.h"
#include "../core/typelist.h"
#include "../core/debug.h"

template<typename T, std::size_t... N>
struct EmptyArray {
  using value_t = T;
  using vector_t = std::vector<value_t>;
  using dimension_t = tags::Dimension<sizeof...(N)>;

  EmptyArray() = default;
  explicit EmptyArray(std::size_t...) {}

  value_t operator()(std::size_t...) const { return static_cast<value_t>(0.0); }
};

template<typename T>
using EmptyArray1D = EmptyArray<T, 1>;

template<typename T>
using EmptyArray2D = EmptyArray<T, 2>;

template<typename T>
using EmptyArray3D = EmptyArray<T, 3>;

template<typename TL>
struct EMDataBase {
  using ex_t = TypeListAt<0, TL>;
  using ey_t = TypeListAt<1, TL>;
  using ez_t = TypeListAt<2, TL>;
  using hx_t = TypeListAt<3, TL>;
  using hy_t = TypeListAt<4, TL>;
  using hz_t = TypeListAt<5, TL>;
  using jx_t = TypeListAt<6, TL>;
  using jy_t = TypeListAt<7, TL>;
  using jz_t = TypeListAt<8, TL>;
  using value_t = typename ex_t::value_t;
  using dimension_t = typename ex_t::dimension_t;

  EMDataBase() = default;

  // explicit EMDataBase(std::size_t nx)
  // : Ex(nx), Ey(nx), Ez(nx),
  //   Hx(nx), Hy(nx), Hz(nx),
  //   Jx(nx), Jy(nx), Jz(nx)
  // {}

  EMDataBase(std::size_t nx, std::size_t ny)
  : Ex(nx, ny), Ey(nx, ny), Ez(nx, ny),
    Hx(nx, ny), Hy(nx, ny), Hz(nx, ny),
    Jx(nx, ny), Jy(nx, ny), Jz(nx, ny)
  {}

  // EMDataBase(std::size_t nx, std::size_t ny, std::size_t nz)
  // : Ex(nx, ny, nz), Ey(nx, ny, nz), Ez(nx, ny, nz),
  //   Hx(nx, ny, nz), Hy(nx, ny, nz), Hz(nx, ny, nz),
  //   Jx(nx, ny, nz), Jy(nx, ny, nz), Jz(nx, ny, nz)
  // {}

  ex_t Ex;
  ey_t Ey;
  ez_t Ez;
  hx_t Hx;
  hy_t Hy;
  hz_t Hz;
  jx_t Jx;
  jy_t Jy;
  jz_t Jz;
};

enum class Derivative{ DX, DY, DZ, NoOp };

template<size_t DIM>
struct IntegratorOffsets {
  std::array<size_t, 2 * DIM> offsets;
};

//====== Curl Operators =======
//=============================
template<Derivative D, bool Forward, typename... IDXS>
struct curl {
  static constexpr auto apply(const auto&, IDXS...) { return 0.0; }
};

template<bool Forward, typename... IDXS>
struct curl<Derivative::DX, Forward, IDXS...> {
  static auto apply(const auto& f, IDXS... idxs) {
    DBG("curl<DX>::apply()", Forward);
    if constexpr (Forward) {
      return f.forward_diff_x(idxs...);
    } else {
      return f.backward_diff_x(idxs...);
    }
  }
};

// template<bool Forward, typename... IDXS>
// struct curl<Derivative::DY, Forward, IDXS...> {
//   static auto apply(const auto& f, IDXS... idxs) { return f.diff_y(idxs...); }
// };

template<bool Forward, typename... IDXS>
struct curl<Derivative::DZ, Forward, IDXS...> {
  static auto apply(const auto& f, IDXS... idxs) {
    DBG("curl<DZ>::apply()", Forward);
    if constexpr (Forward) {
      return f.forward_diff_z(idxs...);
    } else {
      return f.backward_diff_z(idxs...);
    }
  }
};

//=================== Field Functors ========================
//===========================================================
template<typename TL, Derivative ACURL, Derivative BCURL, bool Forward, typename... IDXS>
struct UpdateFunctor {
  using curl1 = curl<ACURL, Forward, IDXS...>;
  using curl2 = curl<BCURL, Forward, IDXS...>;
  using  F = TypeListAt<0, TL>; // Primary Field
  using D1 = TypeListAt<1, TL>; // First Curl Field
  using D2 = TypeListAt<2, TL>; // Second Curl Field
  using J1 = TypeListAt<3, TL>; // Current Density
  using C1 = TypeListAt<4, TL>; // Primary coefficients
  using C2 = TypeListAt<5, TL>; // Curl Coefficients
  using C3 = TypeListAt<6, TL>; // Current Density Coefficients

  using value_t = typename F::value_t;

  static value_t apply(F& f, const D1& d1, const D2& d2, const J1& j, const C1& c_f, const C2& c_d, const C3& c_j, IDXS... idxs) {
    DBG("UpdateFunctor::apply()", ACURL, BCURL);
    const auto    self = c_f(idxs...) * f(idxs...);
    const auto   diff1 = curl1::apply(d1, idxs...);
    const auto   diff2 = curl2::apply(d2, idxs...);
    const auto    diff = c_d(idxs...) * (diff1 - diff2);
    const auto current = c_j(idxs...) * j(idxs...);
    return self + diff - current;
  }
};

// template<typename TL, Derivative ACURL, Derivative BCURL>
// struct FieldIntegrator1D {
//   using dimension_t = tf::tags::Dimension<1>;
//   using offset_t = IntegratorOffsets<1>;
//
//   using  F = TypeListAt<0, TL>; // Primary Field
//   using D1 = TypeListAt<1, TL>; // First Curl Field
//   using D2 = TypeListAt<2, TL>; // Second Curl Field
//   using J1 = TypeListAt<3, TL>; // Current Density
//   using C1 = TypeListAt<4, TL>; // Primary coefficients
//   using C2 = TypeListAt<5, TL>; // Curl Coefficients
//   using C3 = TypeListAt<6, TL>; // Current Density Coefficients
//   using update_func = UpdateFunctor<TL, ACURL, BCURL, std::size_t>;
//   
//   static auto apply(F& f, const D1& d1, const D2& d2, const J1& j1, const C1& c1, const C2& c2, const C3& c3, const offset_t& c4) {
//     for (size_t i = 0; i < f.shape[0]; i++) {
//       update_func::apply(f, d1, d2, j1, c1, c2, c3, c4, i);
//     }
//   }
// };

template<typename TL, Derivative ACURL, Derivative BCURL, bool Forward>
struct FieldIntegrator2D {
  using dimension_t = tf::tags::Dimension<2>;
  using offset_t = IntegratorOffsets<2>;

  using  F = TypeListAt<0, TL>; // Primary Field
  using D1 = TypeListAt<1, TL>; // First Curl Field
  using D2 = TypeListAt<2, TL>; // Second Curl Field
  using J1 = TypeListAt<3, TL>; // Current Density
  using C1 = TypeListAt<4, TL>; // Primary coefficients
  using C2 = TypeListAt<5, TL>; // Curl Coefficients
  using C3 = TypeListAt<6, TL>; // Current Density Coefficients
  using update_func = UpdateFunctor<TL, ACURL, BCURL, Forward, std::size_t, std::size_t>;

  static void apply(F& f, const D1& d1, const D2& d2, const J1& js, const C1& c_f, const C2& c_d, const C3& c_src, const offset_t& o) {
    const auto [x0, x1, y0, y1] = o.offsets;
    DBG("FI2D::apply", ACURL, BCURL);
    for (size_t i = x0; i < f.nx - x1; i++) {
      for (size_t j = y0; j < f.nz - y1; j++) {
        DBG(i, j);
        update_func::apply(f, d1, d2, js, c_f, c_d, c_src, i, j);
      }
    }
  }
};

// template<typename TL, Derivative ACURL, Derivative BCURL>
// struct FieldIntegrator3D {
//   using dimension_t = tf::tags::Dimension<3>;
//   using offset_t = IntegratorOffsets<3>;
//
//   using  F = TypeListAt<0, TL>; // Primary Field
//   using D1 = TypeListAt<1, TL>; // First Curl Field
//   using D2 = TypeListAt<2, TL>; // Second Curl Field
//   using J1 = TypeListAt<3, TL>; // Current Density
//   using C1 = TypeListAt<4, TL>; // Primary coefficients
//   using C2 = TypeListAt<5, TL>; // Curl Coefficients
//   using C3 = TypeListAt<6, TL>; // Current Density Coefficients
//   using update_func = UpdateFunctor<TL, ACURL, BCURL, std::size_t, std::size_t, std::size_t>;
//   
//   static auto apply(F& f, const D1& d1, const D2& d2, const J1& j1, const C1& c1, const C2& c2, const C3& c3, const offset_t& o) {
//     for (size_t i = 0; i < f.shape[0]; i++) {
//       for (size_t j = 0; j < f.shape[1]; j++) {
//         for (size_t k = 0; k < f.shape[2]; k++) {
//           update_func::apply(f, d1, d2, j1, c1, c2, c3, i, j, k);
//         }
//       }
//     }
//   }
// };

struct FieldIntegratorNull {
  using dimension_t = tf::tags::Dimension<0>;
  static constexpr void apply(auto&, auto&, auto&, auto&, auto&, auto&, auto&, auto&) {}
};


template<typename TL, typename EIX, typename EIY, typename EIZ, typename HIX, typename HIY, typename HIZ>
struct Electromagnetics {
  // using value_t = typename TypeListAt<0, TL>::value_t;
  // using dimension_t = typename EIX::dimension_t;
  using empty_t = tf::types::EmptyArray<double, 2>;

  struct emdata_t;

  static constexpr empty_t empty{};

  static constexpr IntegratorOffsets<2> one_offsets{1, 1, 1, 1};
  static constexpr IntegratorOffsets<2> zero_offsets{0, 0, 0, 0};

  static void updateE(emdata_t& emdata) {
    DBG("Electromagnetics::updateE()");
    EIX::apply(emdata.Ex, emdata.Hz, emdata.Hy, emdata.Jx, emdata.c_ex_e, emdata.c_ex_h, emdata.c_jx, zero_offsets);
    EIY::apply(emdata.Ey, emdata.Hx, emdata.Hz, emdata.Jy, emdata.c_ey_e, emdata.c_ey_h, emdata.c_jy, zero_offsets);
    EIZ::apply(emdata.Ez, emdata.Hy, emdata.Hx, emdata.Jz, emdata.c_ez_e, emdata.c_ez_h, emdata.c_jz, one_offsets);
  }

  static void updateH(emdata_t& emdata) {
    DBG("Electromagnetics::updateH()");
    HIX::apply(emdata.Hx, emdata.Ey, emdata.Ez, empty, emdata.c_hx_h, emdata.c_hx_e, empty, zero_offsets);
    HIY::apply(emdata.Hy, emdata.Ez, emdata.Ex, empty, emdata.c_hy_h, emdata.c_hy_e, empty, zero_offsets);
    HIZ::apply(emdata.Hz, emdata.Ex, emdata.Ey, empty, emdata.c_hz_h, emdata.c_hz_e, empty, zero_offsets);
  }

  static void advance(emdata_t& emdata) {
    DBG("Electromagnetics::Advance()");
    updateH(emdata);
    updateE(emdata);
  }
};


template<typename TL, typename EIX, typename EIY, typename EIZ, typename BIX, typename BIY, typename BIZ>
struct Electromagnetics<TL, EIX, EIY, EIZ, BIX, BIY, BIZ>::emdata_t : EMDataBase<TL> {
  using typename EMDataBase<TL>::ex_t;
  using typename EMDataBase<TL>::ey_t;
  using typename EMDataBase<TL>::ez_t;
  using typename EMDataBase<TL>::hx_t;
  using typename EMDataBase<TL>::hy_t;
  using typename EMDataBase<TL>::hz_t;
  using typename EMDataBase<TL>::jx_t;
  using typename EMDataBase<TL>::jy_t;
  using typename EMDataBase<TL>::jz_t;
  using typename EMDataBase<TL>::value_t;
  using typename EMDataBase<TL>::dimension_t;

  emdata_t() = default;

  // explicit emdata_t(std::size_t nx)
  // : EMDataBase<TL>(nx),
  //   c_ex_e(nx), c_ex_h(nx),
  //   c_ey_e(nx), c_ey_h(nx),
  //   c_ez_e(nx), c_ez_h(nx),
  //   c_hx_e(nx), c_hx_h(nx),
  //   c_hy_e(nx), c_hy_h(nx),
  //   c_hz_e(nx), c_hz_h(nx),
  //   c_jx(nx), c_jy(nx), c_jz(nx)
  // {}

  explicit emdata_t(std::size_t nx, std::size_t ny)
  : EMDataBase<TL>(nx, ny),
    c_ex_e(nx, ny), c_ex_h(nx, ny),
    c_ey_e(nx, ny), c_ey_h(nx, ny),
    c_ez_e(nx, ny), c_ez_h(nx, ny),
    c_hx_e(nx, ny), c_hx_h(nx, ny),
    c_hy_e(nx, ny), c_hy_h(nx, ny),
    c_hz_e(nx, ny), c_hz_h(nx, ny),
    c_jx(nx, ny), c_jy(nx, ny), c_jz(nx, ny)
  {}

  // explicit emdata_t(std::size_t nx, std::size_t ny, std::size_t nz)
  // : EMDataBase<TL>(nx, ny, nz),
  //   c_ex_e(nx, ny, nz), c_ex_h(nx, ny, nz),
  //   c_ey_e(nx, ny, nz), c_ey_h(nx, ny, nz),
  //   c_ez_e(nx, ny, nz), c_ez_h(nx, ny, nz),
  //   c_hx_e(nx, ny, nz), c_hx_h(nx, ny, nz),
  //   c_hy_e(nx, ny, nz), c_hy_h(nx, ny, nz),
  //   c_hz_e(nx, ny, nz), c_hz_h(nx, ny, nz),
  //   c_jx(nx, ny, nz), c_jy(nx, ny, nz), c_jz(nx, ny, nz)
  // {}

  ex_t c_ex_e;
  ex_t c_ex_h;
  ey_t c_ey_e;
  ey_t c_ey_h;
  ez_t c_ez_e;
  ez_t c_ez_h;
  hx_t c_hx_e;
  hx_t c_hx_h;
  hy_t c_hy_e;
  hy_t c_hy_h;
  hz_t c_hz_e;
  hz_t c_hz_h;
  jx_t c_jx;
  jy_t c_jy;
  jz_t c_jz;
};

#endif //EM_SOLVER_H
