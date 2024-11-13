//
// Created by cepheid on 6/28/24.
//

#ifndef EM_SOLVER_H
#define EM_SOLVER_H

#include "electromagnetics.param"
// #include "../aydenstuff/array.h"
// #include "../core/debug.h"
// #include "em_data.h"
#include <boundaries.h>

#include "curl_operators.h"
#include "offsets.h"

//=================== Field Functors ========================
//===========================================================
template<Derivative CURL1, Derivative CURL2, bool Forward, typename... IDXS>
struct FieldUpdate {
  using CurlA = curl<CURL1, Forward, IDXS...>;
  using CurlB = curl<CURL2, Forward, IDXS...>;

  static void apply(auto& f, const auto& d1, const auto& d2, const auto& j, const auto& c_f, const auto& c_d, const auto& c_j, IDXS... idxs) {
    const auto    self = c_f(idxs...) * f(idxs...);
    const auto   diff1 = CurlA::apply(d1, idxs...);
    const auto   diff2 = CurlB::apply(d2, idxs...);
    const auto    diff = c_d(idxs...) * (diff1 - diff2);
    const auto current = c_j(idxs...) * j(idxs...);
    f(idxs...) = self + diff - current;
  }
};


/*
 * todo: could these be merged into one Integrator? Make y1/z1 be 1 so the loops all run at least once?
 */
template<typename T, typename UpdateFunctor>
struct FieldIntegrator1D {
  using value_t = typename T::value_t;
  using dimension_t = typename T::dimension_t;
  using array_t = Array1D<value_t>;
  using update_func = UpdateFunctor;

  static auto apply(auto& f, const auto& d1, const auto& d2, const auto& js, const auto& c_f, const auto& c_d, const auto& c_src, const auto& o) {
    for (size_t i = o.x0; i < f.nx() - o.x1; ++i) {
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
    for (size_t i = o.x0; i < f.nx() - o.x1; ++i) {
      for (size_t j = o.y0; j < f.ny() - o.y1; ++j) {
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
#pragma omp parallel for collapse(3)
    for (size_t i = o.x0; i < f.nx() - o.x1; ++i) {
      for (size_t j = o.y0; j < f.ny() - o.y1; ++j) {
        for (size_t k = o.z0; k < f.nz() - o.z1; ++k) {
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


template<typename EXI, typename EYI, typename EZI,
         typename HXI, typename HYI, typename HZI>
struct Electromagnetics {
  using value_t = typename EXI::value_t;
  using dimension_t = typename EXI::dimension_t;
  using empty_t = EmptyArray<value_t, dimension_t::value>;

  static constexpr empty_t empty{};
  static constexpr IntegratorOffsets one_offsets{1, 1, 1, 1, 1, 1};
  static constexpr IntegratorOffsets zero_offsets{0, 0, 0, 0, 0, 0};

  static void updateE(auto& emdata) {
    EXI::apply(emdata.Ex, emdata.Hz, emdata.Hy, emdata.Jx, emdata.Cexe, emdata.Cexh, emdata.Cjx, one_offsets);
    EYI::apply(emdata.Ey, emdata.Hx, emdata.Hz, emdata.Jy, emdata.Ceye, emdata.Ceyh, emdata.Cjy, one_offsets);
    EZI::apply(emdata.Ez, emdata.Hy, emdata.Hx, emdata.Jz, emdata.Ceze, emdata.Cezh, emdata.Cjz, one_offsets);
  }

  static void updateH(auto& emdata) {
    HXI::apply(emdata.Hx, emdata.Ey, emdata.Ez, empty, emdata.Chxh, emdata.Chxe, empty, zero_offsets);
    HYI::apply(emdata.Hy, emdata.Ez, emdata.Ex, empty, emdata.Chyh, emdata.Chye, empty, zero_offsets);
    HZI::apply(emdata.Hz, emdata.Ex, emdata.Ey, empty, emdata.Chzh, emdata.Chze, empty, zero_offsets);
  }

  static void updateE_bcs(auto& emdata, auto& bcdata) {
    // DBG("updateE_bcs()");
    // Boundary<Periodic2D<EMFace::X>>::updateE(bcdata.x0.Ez, emdata.Ez);
    // Boundary<Periodic2D<EMFace::Y>>::updateE(bcdata.y0.Ez, emdata.Ez);

    Boundary<Pml2D<EMFace::X>>::updateE(bcdata.x0.Ez, emdata.Ez, emdata.Hy, emdata.Cezh);
  }

  static void updateH_bcs(auto& emdata, auto& bcdata) {
    // DBG("updateH_bcs()");
    // Boundary<Periodic2D<EMFace::X>>::updateH(bcdata.x0.Hx, emdata.Hx);
    // Boundary<Periodic2D<EMFace::X>>::updateH(bcdata.x0.Hy, emdata.Hy);
    //
    // Boundary<Periodic2D<EMFace::Y>>::updateH(bcdata.y0.Hx, emdata.Hx);
    // Boundary<Periodic2D<EMFace::Y>>::updateH(bcdata.y0.Hy, emdata.Hy);

    // Boundary<Pml2D<EMFace::X>>::updateH(bcdata.x0.Hx, emdata.Hx, emdata.Ez, emdata.Chxe);
    Boundary<Pml2D<EMFace::X>>::updateH(bcdata.x0.Hy, emdata.Hy, emdata.Ez, emdata.Chye);
  }


  static void advance(auto& emdata, auto& bcdata) {
    updateH(emdata);
    updateH_bcs(emdata, bcdata);

    updateE(emdata);
    updateE_bcs(emdata, bcdata);
  }
};

#endif //EM_SOLVER_H
