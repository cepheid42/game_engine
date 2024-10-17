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

enum class Derivative { DX, DY, DZ, NoOp };

struct IntegratorOffsets {
  size_t x0, x1, y0, y1, z0, z1;
};

//====== Curl Operators =======
//=============================
template<Derivative D, bool Forward, typename... IDXS>
struct curl {
  static constexpr auto apply(const auto&, IDXS...) { DBG("curl<NoOp>::apply()"); return 0.0; }
};

template<bool Forward, typename... IDXS>
struct curl<Derivative::DX, Forward, IDXS...> {
  static auto apply(const auto& f, IDXS... idxs) {
    DBG("curl<DX>::apply()");
    if constexpr (Forward) {
      return f.forward_diff_x(idxs...);
    } else {
      return f.backward_diff_x(idxs...);
    }
  }
};

template<bool Forward, typename... IDXS>
struct curl<Derivative::DY, Forward, IDXS...> {
  static auto apply(const auto& f, IDXS... idxs) {
    DBG("curl<DY>::apply()");
    if constexpr (Forward) {
      return f.forward_diff_y(idxs...);
    } else {
      return f.backward_diff_y(idxs...);
    }
  }
};

template<bool Forward, typename... IDXS>
struct curl<Derivative::DZ, Forward, IDXS...> {
  static auto apply(const auto& f, IDXS... idxs) {
    DBG("curl<DZ>::apply()");
    if constexpr (Forward) {
      return f.forward_diff_z(idxs...);
    } else {
      return f.backward_diff_z(idxs...);
    }
  }
};

//=================== Field Functors ========================
//===========================================================
template<Derivative ACURL, Derivative BCURL, bool Forward, typename... IDXS>
struct UpdateFunctor {
  using curl1 = curl<ACURL, Forward, IDXS...>;
  using curl2 = curl<BCURL, Forward, IDXS...>;

  static auto apply(auto& f, const auto& d1, const auto& d2, const auto& j, const auto& c_f, const auto& c_d, const auto& c_j, IDXS... idxs) {
    DBG("UpdateFunctor::apply()");
    const auto    self = c_f(idxs...) * f(idxs...);
    const auto   diff1 = curl1::apply(d1, idxs...);
    const auto   diff2 = curl2::apply(d2, idxs...);
    const auto    diff = c_d(idxs...) * (diff1 - diff2);
    const auto current = c_j(idxs...) * j(idxs...);
    return self + diff - current;
  }
};

template<typename T, Derivative ACURL, Derivative BCURL, bool Forward>
struct FieldIntegrator1D {
  using dimension_t = tf::tags::Dimension<1>;
  using update_func = UpdateFunctor<ACURL, BCURL, Forward, std::size_t>;

  static auto apply(auto& f, const auto& d1, const auto& d2, const auto& js, const auto& c_f, const auto& c_d, const auto& c_src, const auto& o) {
    DBG("FI1D::apply()");
    for (size_t i = o.x0; i < f.nx - o.x1; ++i) {
      update_func::apply(f, d1, d2, js, c_f, c_d, c_src, i);
    }
  }
};

template<typename T, Derivative ACURL, Derivative BCURL, bool Forward>
struct FieldIntegrator2D {
  using value_t = typename T::value_t;
  using dimension_t = typename T::dimension_t;
  using update_func = UpdateFunctor<ACURL, BCURL, Forward, std::size_t, std::size_t>;

  static void apply(auto& f, const auto& d1, const auto& d2, const auto& js, const auto& c_f, const auto& c_d, const auto& c_src, const auto& o) {
    DBG("FI2D::apply()");
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

  static constexpr void apply(auto&, auto&, auto&, auto&, auto&, auto&, auto&, auto&) { DBG("FieldIntegratorNull::apply()"); }
};


template<typename EIX, typename EIY, typename EIZ, typename HIX, typename HIY, typename HIZ>
struct Electromagnetics {
  using value_t = typename EIX::value_t;
  using dimension_t = typename EIX::dimension_t;

  using empty_t = EmptyArray<double, dimension_t::value>;

  static constexpr empty_t empty{};

  static constexpr IntegratorOffsets one_offsets{1, 1, 1, 1, 1, 1};
  static constexpr IntegratorOffsets zero_offsets{0, 0, 0, 0, 0, 0};

  static void updateE(auto& emdata) {
    DBG("Electromagnetics::updateE()");
    EIX::apply(emdata.Ex, emdata.Hz, emdata.Hy, emdata.Jx, emdata.Cexe, emdata.Cexh, emdata.Cjx, one_offsets);
    EIY::apply(emdata.Ey, emdata.Hx, emdata.Hz, emdata.Jy, emdata.Ceye, emdata.Ceyh, emdata.Cjy, one_offsets);
    EIZ::apply(emdata.Ez, emdata.Hy, emdata.Hx, emdata.Jz, emdata.Ceze, emdata.Cezh, emdata.Cjz, one_offsets);
  }

  static void updateH(auto& emdata) {
    DBG("Electromagnetics::updateH()");
    HIX::apply(emdata.Hx, emdata.Ey, emdata.Ez, empty, emdata.Chxh, emdata.Chxe, empty, zero_offsets);
    HIY::apply(emdata.Hy, emdata.Ez, emdata.Ex, empty, emdata.Chyh, emdata.Chye, empty, zero_offsets);
    HIZ::apply(emdata.Hz, emdata.Ex, emdata.Ey, empty, emdata.Chzh, emdata.Chze, empty, zero_offsets);
  }

  static void advance(auto& emdata) {
    DBG("Electromagnetics::Advance()");
    updateH(emdata);
    updateE(emdata);
  }
};

#endif //EM_SOLVER_H
