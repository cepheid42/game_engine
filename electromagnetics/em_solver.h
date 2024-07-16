//
// Created by cepheid on 6/28/24.
//

#ifndef EM_SOLVER_H
#define EM_SOLVER_H

#include <cstdarg>

#include "enums.h"
#include "curls.h"
#include "bcs.h"

using fptype = double;


//=================== Update Functor ========================
//===========================================================
template<typename TL, FieldType FT, Derivative ACURL, Derivative BCURL, typename... IDXS>
struct UpdateFunctor {
  using curl1 = curl<ACURL, FT, IDXS...>;
  using curl2 = curl<BCURL, FT, IDXS...>;
  using F = TypeListAt<0, TL>;
  using D1 = TypeListAt<1, TL>;
  using D2 = TypeListAt<2, TL>;
  using S = TypeListAt<3, TL>;
  using C1 = TypeListAt<4, TL>;
  using C2 = TypeListAt<5, TL>;
  using C3 = TypeListAt<6, TL>;

  static auto apply(F& f, const D1& d1, const D2& d2, const S& src, const C1& c_f, const C2& c_d, const C3& c_src, IDXS... idxs) {
    const auto self = c_f(idxs...) * f(idxs...);
    const auto diff1 = curl1::apply(d1, idxs...);
    const auto diff2 = curl2::apply(d2, idxs...);
    const auto diff = c_d(idxs...) * (diff1 - diff2);
    const auto current = c_src(idxs...) * src(idxs...);
    f(idxs...) = self + diff - current;
  }
};

//=================== Field Integrators ========================
//==============================================================
struct FieldIntegratorNull { static auto apply(...) {} };

template<typename TL, FieldType FT, Derivative ACURL, Derivative BCURL>
struct FieldIntegrator1D {
  using F = TypeListAt<0, TL>;
  using D1 = TypeListAt<1, TL>;
  using D2 = TypeListAt<2, TL>;
  using J1 = TypeListAt<3, TL>;
  using C1 = TypeListAt<4, TL>;
  using C2 = TypeListAt<5, TL>;
  using C3 = TypeListAt<6, TL>;
  using C4 = TypeListAt<7, TL>;
  using update_func = UpdateFunctor<TL, FT, ACURL, BCURL, std::size_t>;
  static auto apply(F& f, const D1& d1, const D2& d2, const J1& j1, const C1& c1, const C2& c2, const C3& c3, const C4& c4) {
    for (size_t i = 0; i < f.shape[0]; i++) {
      update_func::apply(f, d1, d2, j1, c1, c2, c3, c4, i);
    }
  }
};

template<typename TL, FieldType FT, Derivative ACURL, Derivative BCURL>
struct FieldIntegrator2D {
  using F = TypeListAt<0, TL>;
  using D1 = TypeListAt<1, TL>;
  using D2 = TypeListAt<2, TL>;
  using S = TypeListAt<3, TL>;
  using C1 = TypeListAt<4, TL>;
  using C2 = TypeListAt<5, TL>;
  using C3 = TypeListAt<6, TL>;
  using update_func = UpdateFunctor<TL, FT, ACURL, BCURL, std::size_t, std::size_t>;
  static auto apply(F& f, const D1& d1, const D2& d2, const S& src, const C1& c_f, const C2& c_d, const C3& c_src, const IntegratorOffsets<2>& o) {
    const auto [x0, x1, y0, y1] = o.offsets;
    for (size_t i = x0; i < f.shape[0] - x1; i++) {
      for (size_t j = y0; j < f.shape[1] - y1; j++) {
        update_func::apply(f, d1, d2, src, c_f, c_d, c_src, i, j);
      }
    }
  }
};

template<typename TL, FieldType FT, Derivative ACURL, Derivative BCURL>
struct FieldIntegrator3D {
  using F = TypeListAt<0, TL>;
  using D1 = TypeListAt<1, TL>;
  using D2 = TypeListAt<2, TL>;
  using S = TypeListAt<3, TL>;
  using C1 = TypeListAt<4, TL>;
  using C2 = TypeListAt<5, TL>;
  using C3 = TypeListAt<6, TL>;
  using update_func = UpdateFunctor<TL, FT, ACURL, BCURL, std::size_t, std::size_t, std::size_t>;
  static auto apply(F& f, const D1& d1, const D2& d2, const S& src, const C1& c_f, const C2& c_d, const C3& c_src, const IntegratorOffsets<3>& o) {
    const auto [x0, x1, y0, y1, z0, z1] = o.offsets;
    for (size_t i = x0; i < f.shape[0] - x1; i++) {
      for (size_t j = y0; j < f.shape[1] - y1; j++) {
        for (size_t k = z0; k < f.shape[2] - z1; k++) {
          update_func::apply(f, d1, d2, src, c_f, c_d, c_src, i, j, k);
        }
      }
    }
  }
};

//=================== Electromagnetics Class ========================
//===================================================================
template<typename EIX, typename EIY, typename EIZ, typename HIX, typename HIY, typename HIZ>
struct Electromagnetics {
  static void updateE(auto& fields, const auto& coeffs) {
    constexpr IntegratorOffsets<2> ex_offsets = {0u, 0u, 0u, 0u};
    constexpr IntegratorOffsets<2> ey_offsets = {0u, 0u, 0u, 0u};
    constexpr IntegratorOffsets<2> ez_offsets = {1u, 1u, 1u, 1u};

    // constexpr IntegratorOffsets<3> ex_offsets = {0u, 0u, 1u, 1u, 1u, 1u};
    // constexpr IntegratorOffsets<3> ey_offsets = {1u, 1u, 0u, 0u, 1u, 1u};
    // constexpr IntegratorOffsets<3> ez_offsets = {1u, 1u, 1u, 1u, 0u, 0u};

    EIX::apply(fields.Ex, fields.Hz, fields.Hy, fields.Jx, coeffs.cexe, coeffs.cexh, coeffs.cjx, ex_offsets);
    EIY::apply(fields.Ey, fields.Hx, fields.Hz, fields.Jy, coeffs.ceye, coeffs.ceyh, coeffs.cjy, ey_offsets);
    EIZ::apply(fields.Ez, fields.Hy, fields.Hx, fields.Jz, coeffs.ceze, coeffs.cezh, coeffs.cjz, ez_offsets);
  }

  static void updateB(auto& fields, const auto& coeffs) {
    constexpr IntegratorOffsets<2> hx_offsets = {0u, 0u, 0u, 0u};
    constexpr IntegratorOffsets<2> hy_offsets = {0u, 0u, 0u, 0u};
    constexpr IntegratorOffsets<2> hz_offsets = {0u, 0u, 0u, 0u};
    // constexpr IntegratorOffsets<3> hx_offsets = {0u, 0u, 0u, 0u, 0u, 0u};
    // constexpr IntegratorOffsets<3> hy_offsets = {0u, 0u, 0u, 0u, 0u, 0u};
    // constexpr IntegratorOffsets<3> hz_offsets = {0u, 0u, 0u, 0u, 0u, 0u};

    HIX::apply(fields.Hx, fields.Ey, fields.Ez, EMPTYARRAY{}, coeffs.chxh, coeffs.chxe, EMPTYARRAY{}, hx_offsets);
    HIY::apply(fields.Hy, fields.Ez, fields.Ex, EMPTYARRAY{}, coeffs.chyh, coeffs.chye, EMPTYARRAY{}, hy_offsets);
    HIZ::apply(fields.Hz, fields.Ex, fields.Ey, EMPTYARRAY{}, coeffs.chzh, coeffs.chze, EMPTYARRAY{}, hz_offsets);
  }
};


#endif //EM_SOLVER_H
