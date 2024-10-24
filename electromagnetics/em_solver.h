//
// Created by cepheid on 6/28/24.
//

#ifndef EM_SOLVER_H
#define EM_SOLVER_H

#include <cassert>
#include <array>


#include "../aydenstuff/array.h"
#include "../core/typelist.h"
#include "../core/debug.h"
#include "em_data.h"
#include "curl_operators.h"
#include "bc_data.h"
#include "boundaries.h"



//=================== Field Functors ========================
//===========================================================
template<Derivative ACURL, Derivative BCURL, bool Forward, typename... IDXS>
struct FieldUpdate {
  using curl1 = curl<ACURL, Forward, IDXS...>;
  using curl2 = curl<BCURL, Forward, IDXS...>;

  static void apply(auto& f, const auto& d1, const auto& d2, const auto& j, const auto& c_f, const auto& c_d, const auto& c_j, IDXS... idxs) {
    // DBG("UpdateFunctor::apply()");
    const auto    self = c_f(idxs...) * f(idxs...);
    const auto   diff1 = curl1::apply(d1, idxs...);
    const auto   diff2 = curl2::apply(d2, idxs...);
    const auto    diff = c_d(idxs...) * (diff1 - diff2);
    const auto current = c_j(idxs...) * j(idxs...);
    // (..., DBG(idxs));
    f(idxs...) = self + diff - current;
  }
};

template<typename T, typename UpdateFunctor>
struct FieldIntegrator1D {
  using value_t = typename T::value_t;
  using dimension_t = typename T::dimension_t;
  using array_t = Array1D<value_t>;
  using update_func = UpdateFunctor;

  static auto apply(auto& f, const auto& d1, const auto& d2, const auto& js, const auto& c_f, const auto& c_d, const auto& c_src, const auto& o) {
    // DBG("FI1D::apply()", o.x0, o.x1, f.nx - o.x1);
    for (size_t i = o.x0; i < f.nx - o.x1; ++i) {
      update_func::apply(f, d1, d2, js, c_f, c_d, c_src, i);
    }
  }
};

template<typename T, typename UpdateFunctor>
struct FieldIntegrator2D {
  using value_t = typename T::value_t;
  using dimension_t = typename T::dimension_t;
  using array_t = Array2D<value_t>;
  using update_func = UpdateFunctor;

  static void apply(auto& f, const auto& d1, const auto& d2, const auto& js, const auto& c_f, const auto& c_d, const auto& c_src, const auto& o) {
    // DBG("FI2D::apply()");
    for (size_t i = o.x0; i < f.nx - o.x1; ++i) {
      for (size_t j = o.y0; j < f.nz - o.y1; ++j) {
        // DBG(i, j);
        update_func::apply(f, d1, d2, js, c_f, c_d, c_src, i, j);
      }
    }
  }
};

// template<typename TL, Derivative ACURL, Derivative BCURL>
// struct FieldIntegrator3D {
//   using dimension_t = tf::tags::Dimension<3>;
//   using offset_t = IntegratorOffsets<3>;
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

template<typename T>
struct FieldIntegratorNull {
  using value_t = typename T::value_t;
  using dimension_t = typename T::dimension_t;
  using array_t = EmptyArray<value_t, dimension_t::value>;

  static constexpr void apply(const auto&, const auto&, const auto&, const auto&, const auto&, const auto&, const auto&, const auto&) {}
};

struct IntegratorOffsets {
  size_t x0, x1, y0, y1, z0, z1;
};

template<typename EIX, typename EIY, typename EIZ, typename HIX, typename HIY, typename HIZ, typename LOBC, typename HIBC>
struct Electromagnetics {
  using value_t = typename EIX::value_t;
  using dimension_t = typename EIX::dimension_t;
  using empty_t = EmptyArray<value_t, dimension_t::value>;

  static constexpr empty_t empty{};
  static constexpr IntegratorOffsets Eoffsets{1, 1, 1, 1, 1, 1};
  static constexpr IntegratorOffsets Hoffsets{0, 0, 0, 0, 0, 0};

  static constexpr IntegratorOffsets ExPeriodic{0, 0, 0, 0, 0, 0};
  static constexpr IntegratorOffsets EyPeriodic{0, 0, 0, 0, 0, 0};
  static constexpr IntegratorOffsets EzPeriodic{0, nHalo, 0, 0, 0, 0};
  static constexpr IntegratorOffsets HxPeriodic{0, 0, 0, 0, 0, 0};
  static constexpr IntegratorOffsets HyPeriodic{0, nHalo, 0, 0, 0, 0};
  static constexpr IntegratorOffsets HzPeriodic{0, 0, 0, 0, 0, 0};

  using ExLo = TypeListAt<0, LOBC>;
  using EyLo = TypeListAt<1, LOBC>;
  using EzLo = TypeListAt<2, LOBC>;
  using HxLo = TypeListAt<3, LOBC>;
  using HyLo = TypeListAt<4, LOBC>;
  using HzLo = TypeListAt<5, LOBC>;
  
  using ExHi = TypeListAt<0, HIBC>;
  using EyHi = TypeListAt<1, HIBC>;
  using EzHi = TypeListAt<2, HIBC>;
  using HxHi = TypeListAt<3, HIBC>;
  using HyHi = TypeListAt<4, HIBC>;
  using HzHi = TypeListAt<5, HIBC>;

  static void updateELoBC(auto& emdata, auto& bcdata) {
    DBG("Electromagnetics::updateELoBC()");
    ExLo::apply(emdata.Ex, emdata.Hz, emdata.Hy, emdata.Cexh, bcdata.psiEx, bcdata.bEx, bcdata.cEx, ExPeriodic);
    EyLo::apply(emdata.Ey, emdata.Hx, emdata.Hz, emdata.Ceyh, bcdata.psiEy, bcdata.bEy, bcdata.cEy, EyPeriodic);
    EzLo::apply(emdata.Ez, emdata.Hy, emdata.Hx, emdata.Cezh, bcdata.psiEz, bcdata.bEz, bcdata.cEz, EzPeriodic);
  }

  static void updateHLoBC(auto& emdata, auto& bcdata) {
    DBG("Electromagnetics::updateHLoBC()");
    HxLo::apply(emdata.Hx, emdata.Ey, emdata.Ez, emdata.Chxe, bcdata.psiHx, bcdata.bHx, bcdata.cHx, HxPeriodic);
    HyLo::apply(emdata.Hy, emdata.Ez, emdata.Ex, emdata.Chye, bcdata.psiHy, bcdata.bHy, bcdata.cHy, HyPeriodic);
    HzLo::apply(emdata.Hz, emdata.Ex, emdata.Ey, emdata.Chze, bcdata.psiHz, bcdata.bHz, bcdata.cHz, HzPeriodic);
  }

  static void updateE(auto& emdata, auto& bcdata) {
    // DBG("Electromagnetics::updateE()");
    EIX::apply(emdata.Ex, emdata.Hz, emdata.Hy, emdata.Jx, emdata.Cexe, emdata.Cexh, emdata.Cjx, Eoffsets);
    EIY::apply(emdata.Ey, emdata.Hx, emdata.Hz, emdata.Jy, emdata.Ceye, emdata.Ceyh, emdata.Cjy, Eoffsets);
    EIZ::apply(emdata.Ez, emdata.Hy, emdata.Hx, emdata.Jz, emdata.Ceze, emdata.Cezh, emdata.Cjz, Eoffsets);

    updateELoBC(emdata, bcdata);
    // ExLo::apply(emdata.Ex, emdata.Hz, emdata.Hy, emdata.Cexh, bcdata.psiEx, bcdata.bEx, bcdata.cEx, ExPeriodic);
    // EyLo::apply(emdata.Ey, emdata.Hx, emdata.Hz, emdata.Ceyh, bcdata.psiEy, bcdata.bEy, bcdata.cEy, EyPeriodic);
    // EzLo::apply(emdata.Ez, emdata.Hy, emdata.Hx, emdata.Cezh, bcdata.psiEz, bcdata.bEz, bcdata.cEz, EzPeriodic);

    // ExHi::apply(emdata.Ex, emdata.Hz, emdata.Hy, emdata.Cexh, bcdata.psiEx, bcdata.bEx, bcdata.cEx, ExPeriodic);
    // EyHi::apply(emdata.Ey, emdata.Hx, emdata.Hz, emdata.Ceyh, bcdata.psiEy, bcdata.bEy, bcdata.cEy, EyPeriodic);
    // EzHi::apply(emdata.Ez, emdata.Hy, emdata.Hx, emdata.Cezh, bcdata.psiEz, bcdata.bEz, bcdata.cEz, EzPeriodic);
  }

  static void updateH(auto& emdata, auto& bcdata) {
    // DBG("Electromagnetics::updateH()");
    HIX::apply(emdata.Hx, emdata.Ey, emdata.Ez, empty, emdata.Chxh, emdata.Chxe, empty, Hoffsets);
    HIY::apply(emdata.Hy, emdata.Ez, emdata.Ex, empty, emdata.Chyh, emdata.Chye, empty, Hoffsets);
    HIZ::apply(emdata.Hz, emdata.Ex, emdata.Ey, empty, emdata.Chzh, emdata.Chze, empty, Hoffsets);

    updateHLoBC(emdata, bcdata);
    // HxLo::apply(emdata.Hx, emdata.Ey, emdata.Ez, emdata.Chxe, bcdata.psiHx, bcdata.bHx, bcdata.cHx, HxPeriodic);
    // HyLo::apply(emdata.Hy, emdata.Ez, emdata.Ex, emdata.Chye, bcdata.psiHy, bcdata.bHy, bcdata.cHy, HyPeriodic);
    // HzLo::apply(emdata.Hz, emdata.Ex, emdata.Ey, emdata.Chze, bcdata.psiHz, bcdata.bHz, bcdata.cHz, HzPeriodic);

    // HxHi::apply(emdata.Hx, emdata.Ey, emdata.Ez, emdata.Chxe, bcdata.psiHx, bcdata.bHx, bcdata.cHx, HxPeriodic);
    // HyHi::apply(emdata.Hy, emdata.Ez, emdata.Ex, emdata.Chye, bcdata.psiHy, bcdata.bHy, bcdata.cHy, HyPeriodic);
    // HzHi::apply(emdata.Hz, emdata.Ex, emdata.Ey, emdata.Chze, bcdata.psiHz, bcdata.bHz, bcdata.cHz, HzPeriodic);
  }


  static void advance(auto& emdata, auto& bcdata) {
    // DBG("Electromagnetics::Advance()");
    updateH(emdata, bcdata);
    updateE(emdata, bcdata);
  }
};

#endif //EM_SOLVER_H
