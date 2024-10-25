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
#include "offsets.h"
// #include "bc_data.h"
// #include "boundaries.h"



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


// template<EMFace F, bool=false>
// struct PeriodicOffsets {
//
//   explicit PeriodicOffsets(const auto& f) requires (F == EMFace::X)
//   : offsets{0, nHalo, 0, f.ny, 0, f.nz}
//   {}
//
//   explicit PeriodicOffsets(const auto& f) requires (F == EMFace::Y)
//   : offsets{0, f.nx, 0, nHalo, 0, f.nz}
//   {}
//
//   explicit PeriodicOffsets(const auto& f) requires (F == EMFace::Z)
//   : offsets{0, f.nx, 0, f.ny, 0, nHalo}
//   {}
//
//   explicit PeriodicOffsets(auto&) {}
//
//   IntegratorOffsets offsets{};
// };
//
// template<EMFace F, bool HI>
// struct PMLOffsets {
//
//   explicit PMLOffsets(const auto& f) requires (F == EMFace::X and !HI)
//   : offsets{0, dPML, 0, f.ny, 0, f.nz}
//   {}
//
//   explicit PMLOffsets(const auto& f) requires (F == EMFace::Y and !HI)
//   : offsets{0, f.nx, 0, nHalo, 0, f.nz}
//   {}
//
//   explicit PMLOffsets(const auto& f) requires (F == EMFace::Z and !HI)
//   : offsets{0, f.nx, 0, f.ny, 0, nHalo}
//   {}
//
//   explicit PMLOffsets(const auto& f) requires (F == EMFace::X and HI)
//   : offsets{f.nx - dPML, f.nx, 0, f.ny, 0, f.nz}
//   {}
//
//   explicit PMLOffsets(const auto& f) requires (F == EMFace::Y and HI)
//   : offsets{0, f.nx, 0, nHalo, 0, f.nz}
//   {}
//
//   explicit PMLOffsets(const auto& f) requires (F == EMFace::Z and HI)
//   : offsets{0, f.nx, 0, f.ny, 0, nHalo}
//   {}
//
//   explicit PMLOffsets(auto&) {}
//
//   IntegratorOffsets offsets{};
// };


template<typename EIX, typename EIY, typename EIZ, typename HIX, typename HIY, typename HIZ, typename LOBC>//, typename HIBC>
struct Electromagnetics {
  using value_t = typename EIX::value_t;
  using dimension_t = typename EIX::dimension_t;
  using empty_t = EmptyArray<value_t, dimension_t::value>;

  static constexpr empty_t empty{};
  static constexpr IntegratorOffsets one_offsets{1, 1, 1, 1, 1, 1};
  static constexpr IntegratorOffsets zero_offsets{0, 0, 0, 0, 0, 0};

  using XLo = TypeListAt<0, LOBC>;
  using YLo = TypeListAt<1, LOBC>;
  using ZLo = TypeListAt<2, LOBC>;

  
  // using ExHi = TypeListAt<0, HIBC>;
  // using EyHi = TypeListAt<1, HIBC>;
  // using EzHi = TypeListAt<2, HIBC>;
  // using HxHi = TypeListAt<3, HIBC>;
  // using HyHi = TypeListAt<4, HIBC>;
  // using HzHi = TypeListAt<5, HIBC>;

  static void updateELoBC(auto& emdata, auto& bcdata) {
    DBG("Electromagnetics::updateELoBC()");
    // const auto Ex_offsets =
    // XLo::apply(emdata.Ex, emdata.Hz, emdata.Hy, emdata.Cexh, bcdata.psiEx, bcdata.bEx, bcdata.cEx, {0, 0, 0, 0, 0, 0});
    // EyLo::apply(emdata.Ey, emdata.Hx, emdata.Hz, emdata.Ceyh, bcdata.psiEy, bcdata.bEy, bcdata.cEy, EyPeriodic);
    // EzLo::apply(emdata.Ez, emdata.Hy, emdata.Hx, emdata.Cezh, bcdata.psiEz, bcdata.bEz, bcdata.cEz, EzPeriodic);
  }

  static void updateHLoBC(auto& emdata, auto& bcdata) {
    DBG("Electromagnetics::updateHLoBC()");
    // HxLo::apply(emdata.Hx, emdata.Ey, emdata.Ez, emdata.Chxe, bcdata.psiHx, bcdata.bHx, bcdata.cHx, HxPeriodic);
    // HyLo::apply(emdata.Hy, emdata.Ez, emdata.Ex, emdata.Chye, bcdata.psiHy, bcdata.bHy, bcdata.cHy, HyPeriodic);
    // HzLo::apply(emdata.Hz, emdata.Ex, emdata.Ey, emdata.Chze, bcdata.psiHz, bcdata.bHz, bcdata.cHz, HzPeriodic);
  }

  static void updateE(auto& emdata, auto& bcdata) {
    // DBG("Electromagnetics::updateE()");
    EIX::apply(emdata.Ex, emdata.Hz, emdata.Hy, emdata.Jx, emdata.Cexe, emdata.Cexh, emdata.Cjx, one_offsets);
    EIY::apply(emdata.Ey, emdata.Hx, emdata.Hz, emdata.Jy, emdata.Ceye, emdata.Ceyh, emdata.Cjy, one_offsets);
    EIZ::apply(emdata.Ez, emdata.Hy, emdata.Hx, emdata.Jz, emdata.Ceze, emdata.Cezh, emdata.Cjz, one_offsets);

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
    HIX::apply(emdata.Hx, emdata.Ey, emdata.Ez, empty, emdata.Chxh, emdata.Chxe, empty, zero_offsets);
    HIY::apply(emdata.Hy, emdata.Ez, emdata.Ex, empty, emdata.Chyh, emdata.Chye, empty, zero_offsets);
    HIZ::apply(emdata.Hz, emdata.Ex, emdata.Ey, empty, emdata.Chzh, emdata.Chze, empty, zero_offsets);

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
    // DBG(dbg::type<LOBC>());
    DBG(dbg::type<XLo>(), dbg::type<YLo>(), dbg::type<ZLo>());
    // updateH(emdata, bcdata);
    // updateE(emdata, bcdata);
  }
};

#endif //EM_SOLVER_H
