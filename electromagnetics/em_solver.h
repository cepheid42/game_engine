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

  static constexpr auto dx = 1.0 / 99.0;

  static void apply(auto& f, const auto& d1, const auto& d2, const auto& j, const auto& c_f, const auto& c_d, const auto& c_j, IDXS... idxs) {
    const auto    self = c_f(idxs...) * f(idxs...);
    const auto   diff1 = CurlA::apply(d1, idxs...);
    const auto   diff2 = CurlB::apply(d2, idxs...);
    const auto    diff = (c_d(idxs...) / dx) * (diff1 - diff2);
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
  using offset_t = IntegratorOffsets;

  static auto apply(auto& f, const auto& d1, const auto& d2, const auto& js, const auto& c_f, const auto& c_d, const auto& c_src, const offset_t& o) {
#pragma omp parallel for num_threads(NTHREADS)
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
  using offset_t = IntegratorOffsets;

  static void apply(auto& f, const auto& d1, const auto& d2, const auto& js, const auto& c_f, const auto& c_d, const auto& c_src, const offset_t& o) {
#pragma omp parallel for collapse(2) num_threads(NTHREADS)
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
  using offset_t = IntegratorOffsets;

  static auto apply(auto& f, const auto& d1, const auto& d2, const auto& js, const auto& c_f, const auto& c_d, const auto& c_src, const offset_t& o) {
#pragma omp parallel for collapse(3) num_threads(NTHREADS)
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
  using offset_t = IntegratorOffsets;

  static constexpr void apply(const auto&, const auto&, const auto&, const auto&, const auto&, const auto&, const auto&, const offset_t&) {}
};


template<typename EXI, typename EYI, typename EZI,
         typename HXI, typename HYI, typename HZI,
         typename BCX0, typename BCX1,
         typename BCY0, typename BCY1,
         typename BCZ0, typename BCZ1>
struct Electromagnetics {
  using value_t = typename EXI::value_t;
  using dimension_t = typename EXI::dimension_t;
  using empty_t = EmptyArray<value_t, dimension_t::value>;

  static constexpr empty_t empty{};
  // static constexpr IntegratorOffsets zero_offsets{0, 0, 0, 0, 0, 0};

  static void updateE(auto& emdata) {
    EXI::apply(emdata.Ex, emdata.Hz, emdata.Hy, emdata.Jx, emdata.Cexe, emdata.Cexh, emdata.Cjx, {0, 0, 1, 1, 1, 1});
    EYI::apply(emdata.Ey, emdata.Hx, emdata.Hz, emdata.Jy, emdata.Ceye, emdata.Ceyh, emdata.Cjy, {1, 1, 0, 0, 1, 1});
    EZI::apply(emdata.Ez, emdata.Hy, emdata.Hx, emdata.Jz, emdata.Ceze, emdata.Cezh, emdata.Cjz, {1, 1, 1, 1, 0, 0});
  }

  static void updateH(auto& emdata) {
    HXI::apply(emdata.Hx, emdata.Ey, emdata.Ez, empty, emdata.Chxh, emdata.Chxe, empty, {0, 0, 0, 0, 0, 0});
    HYI::apply(emdata.Hy, emdata.Ez, emdata.Ex, empty, emdata.Chyh, emdata.Chye, empty, {0, 0, 0, 0, 0, 0});
    HZI::apply(emdata.Hz, emdata.Ex, emdata.Ey, empty, emdata.Chzh, emdata.Chze, empty, {0, 0, 0, 0, 0, 0});
  }

  static void updateE_bcs(auto& emdata, auto& bcdata) {
    BCX0::Ex::updateE(bcdata.x0.Ex, emdata.Ex, empty, empty);
    BCX0::Ey::updateE(bcdata.x0.Ey, emdata.Ey, emdata.Hz, emdata.Ceyh);
    BCX0::Ez::updateE(bcdata.x0.Ez, emdata.Ez, emdata.Hy, emdata.Cezh);

    BCX1::Ex::updateE(bcdata.x1.Ex, emdata.Ex, empty, empty);
    BCX1::Ey::updateE(bcdata.x1.Ey, emdata.Ey, emdata.Hz, emdata.Ceyh);
    BCX1::Ez::updateE(bcdata.x1.Ez, emdata.Ez, emdata.Hy, emdata.Cezh);

    BCY0::Ex::updateE(bcdata.y0.Ex, emdata.Ex, emdata.Hz, emdata.Cexh);
    BCY0::Ey::updateE(bcdata.y0.Ey, emdata.Ey, empty, empty);
    BCY0::Ez::updateE(bcdata.y0.Ez, emdata.Ez, emdata.Hx, emdata.Cezh);

    BCY1::Ex::updateE(bcdata.y1.Ex, emdata.Ex, emdata.Hz, emdata.Cexh);
    BCY1::Ey::updateE(bcdata.y1.Ey, emdata.Ey, empty, empty);
    BCY1::Ez::updateE(bcdata.y1.Ez, emdata.Ez, emdata.Hx, emdata.Cezh);

    BCZ0::Ex::updateE(bcdata.z0.Ex, emdata.Ex, emdata.Hy, emdata.Cexh);
    BCZ0::Ey::updateE(bcdata.z0.Ey, emdata.Ey, emdata.Hx, emdata.Ceyh);
    BCZ0::Ez::updateE(bcdata.z0.Ez, emdata.Ez, empty, empty);

    BCZ1::Ex::updateE(bcdata.z1.Ex, emdata.Ex, emdata.Hy, emdata.Cexh);
    BCZ1::Ey::updateE(bcdata.z1.Ey, emdata.Ey, emdata.Hx, emdata.Ceyh);
    BCZ1::Ez::updateE(bcdata.z1.Ez, emdata.Ez, empty, empty);
  }

  static void updateH_bcs(auto& emdata, auto& bcdata) {
    BCX0::Hx::updateH(bcdata.x0.Hx, emdata.Hx, empty, empty);
    BCX0::Hy::updateH(bcdata.x0.Hy, emdata.Hy, emdata.Ez, emdata.Chye);
    BCX0::Hz::updateH(bcdata.x0.Hz, emdata.Hz, emdata.Ey, emdata.Chze);

    BCX1::Hx::updateH(bcdata.x1.Hx, emdata.Hx, empty, empty);
    BCX1::Hy::updateH(bcdata.x1.Hy, emdata.Hy, emdata.Ez, emdata.Chye);
    BCX1::Hz::updateH(bcdata.x1.Hz, emdata.Hz, emdata.Ey, emdata.Chze);

    BCY0::Hx::updateH(bcdata.y0.Hx, emdata.Hx, emdata.Ez, emdata.Chxe);
    BCY0::Hy::updateH(bcdata.y0.Hy, emdata.Hy, empty, empty);
    BCY0::Hz::updateH(bcdata.y0.Hz, emdata.Hz, emdata.Ex, emdata.Chze);

    BCY1::Hx::updateH(bcdata.y1.Hx, emdata.Hx, emdata.Ez, emdata.Chxe);
    BCY1::Hy::updateH(bcdata.y1.Hy, emdata.Hy, empty, empty);
    BCY1::Hz::updateH(bcdata.y1.Hz, emdata.Hz, emdata.Ex, emdata.Chze);

    BCZ0::Hx::updateH(bcdata.z0.Hx, emdata.Hx, emdata.Ey, emdata.Chxe);
    BCZ0::Hy::updateH(bcdata.z0.Hy, emdata.Hy, emdata.Ex, emdata.Chye);
    BCZ0::Hz::updateH(bcdata.z0.Hz, emdata.Hz, empty, empty);

    BCZ1::Hx::updateH(bcdata.z1.Hx, emdata.Hx, emdata.Ey, emdata.Chxe);
    BCZ1::Hy::updateH(bcdata.z1.Hy, emdata.Hy, emdata.Ex, emdata.Chye);
    BCZ1::Hz::updateH(bcdata.z1.Hz, emdata.Hz, empty, empty);
  }


  static void advance(auto& emdata, auto& bcdata) {
    updateH(emdata);
    updateH_bcs(emdata, bcdata);

    updateE(emdata);
    updateE_bcs(emdata, bcdata);
  }
};

#endif //EM_SOLVER_H
