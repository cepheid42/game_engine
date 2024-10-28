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

template<typename T, typename UpdateFunctor>
struct FieldIntegrator3D {
  using value_t = typename T::value_t;
  using dimension_t = typename T::dimension_t;
  using array_t = Array3D<value_t>;
  using update_func = UpdateFunctor;

  static auto apply(auto& f, const auto& d1, const auto& d2, const auto& js, const auto& c_f, const auto& c_d, const auto& c_src, const auto& o) {
    for (size_t i = o.x0; i < f.shape[0] - o.x1; ++i) {
      for (size_t j = o.y0; j < f.shape[1] - o.y1; ++j) {
        for (size_t k = o.z0; k < f.shape[2] - o.z1; ++k) {
          update_func::apply(f, d1, d2, js, c_f, c_d, c_src, i, j, k);
        }
      }
    }
  }
};

template<typename T>
struct FieldIntegratorNull {
  using value_t = typename T::value_t;
  using dimension_t = typename T::dimension_t;
  using array_t = EmptyArray<value_t, dimension_t::value>;

  static constexpr void apply(const auto&, const auto&, const auto&, const auto&, const auto&, const auto&, const auto&, const auto&) {}
};






template<typename EIX, typename EIY, typename EIZ, typename HIX, typename HIY, typename HIZ, typename X0BC>//, typename Y0BC, typename Z0BC>
struct Electromagnetics {
  using value_t = typename EIX::value_t;
  using dimension_t = typename EIX::dimension_t;
  using empty_t = EmptyArray<value_t, dimension_t::value>;

  static constexpr empty_t empty{};
  static constexpr IntegratorOffsets one_offsets{1, 1, 1, 1, 1, 1};
  static constexpr IntegratorOffsets zero_offsets{0, 0, 0, 0, 0, 0};

  static void updateEBC(auto& emdata, auto& bcdata) {
    DBG("Electromagnetics::updateELoBC()");
    X0BC::applyE(emdata, bcdata);
    // Y0BC::applyE(emdata, bcdata);
    // Z0BC::applyE(emdata, bcdata);
  }

  static void updateHBC(auto& emdata, auto& bcdata) {
    DBG("Electromagnetics::updateHLoBC()");
    X0BC::applyH(emdata, bcdata);
    // Y0BC::applyH(emdata, bcdata);
    // Z0BC::applyH(emdata, bcdata);
  }

  static void updateE(auto& emdata, auto& bcdata) {
    DBG("Electromagnetics::updateE()");
    EIX::apply(emdata.Ex, emdata.Hz, emdata.Hy, emdata.Jx, emdata.Cexe, emdata.Cexh, emdata.Cjx, one_offsets);
    EIY::apply(emdata.Ey, emdata.Hx, emdata.Hz, emdata.Jy, emdata.Ceye, emdata.Ceyh, emdata.Cjy, one_offsets);
    EIZ::apply(emdata.Ez, emdata.Hy, emdata.Hx, emdata.Jz, emdata.Ceze, emdata.Cezh, emdata.Cjz, one_offsets);

    updateEBC(emdata, bcdata);
  }

  static void updateH(auto& emdata, auto& bcdata) {
    DBG("Electromagnetics::updateH()");
    HIX::apply(emdata.Hx, emdata.Ey, emdata.Ez, empty, emdata.Chxh, emdata.Chxe, empty, zero_offsets);
    HIY::apply(emdata.Hy, emdata.Ez, emdata.Ex, empty, emdata.Chyh, emdata.Chye, empty, zero_offsets);
    HIZ::apply(emdata.Hz, emdata.Ex, emdata.Ey, empty, emdata.Chzh, emdata.Chze, empty, zero_offsets);

    updateHBC(emdata, bcdata);
  }


  static void advance(auto& emdata, auto& bcdata) {
    DBG("Electromagnetics::Advance()");
    updateH(emdata, bcdata);
    updateE(emdata, bcdata);
  }
};

#endif //EM_SOLVER_H
