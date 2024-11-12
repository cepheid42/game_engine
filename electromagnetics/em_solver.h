//
// Created by cepheid on 6/28/24.
//

#ifndef EM_SOLVER_H
#define EM_SOLVER_H

#include "electromagnetics.param"
// #include "../aydenstuff/array.h"
// #include "../core/debug.h"
// #include "em_data.h"
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
         typename HXI, typename HYI, typename HZI,
         typename X0BC, typename X1BC,
         typename Y0BC, typename Y1BC,
         typename Z0BC, typename Z1BC>
struct Electromagnetics {
  using value_t = typename EXI::value_t;
  using dimension_t = typename EXI::dimension_t;
  using empty_t = EmptyArray<value_t, dimension_t::value>;

  static constexpr empty_t empty{};
  static constexpr IntegratorOffsets one_offsets{1, 1, 1, 1, 1, 1};
  static constexpr IntegratorOffsets zero_offsets{0, 0, 0, 0, 0, 0};

  using ExX0BC = TypeListAt<0, X0BC>;
  using EyX0BC = TypeListAt<1, X0BC>;
  using EzX0BC = TypeListAt<2, X0BC>;
  using HxX0BC = TypeListAt<3, X0BC>;
  using HyX0BC = TypeListAt<4, X0BC>;
  using HzX0BC = TypeListAt<5, X0BC>;

  using ExX1BC = TypeListAt<0, X1BC>;
  using EyX1BC = TypeListAt<1, X1BC>;
  using EzX1BC = TypeListAt<2, X1BC>;
  using HxX1BC = TypeListAt<3, X1BC>;
  using HyX1BC = TypeListAt<4, X1BC>;
  using HzX1BC = TypeListAt<5, X1BC>;

  using ExY0BC = TypeListAt<0, Y0BC>;
  using EyY0BC = TypeListAt<1, Y0BC>;
  using EzY0BC = TypeListAt<2, Y0BC>;
  using HxY0BC = TypeListAt<3, Y0BC>;
  using HyY0BC = TypeListAt<4, Y0BC>;
  using HzY0BC = TypeListAt<5, Y0BC>;

  using ExY1BC = TypeListAt<0, Y1BC>;
  using EyY1BC = TypeListAt<1, Y1BC>;
  using EzY1BC = TypeListAt<2, Y1BC>;
  using HxY1BC = TypeListAt<3, Y1BC>;
  using HyY1BC = TypeListAt<4, Y1BC>;
  using HzY1BC = TypeListAt<5, Y1BC>;

  using ExZ0BC = TypeListAt<0, Z0BC>;
  using EyZ0BC = TypeListAt<1, Z0BC>;
  using EzZ0BC = TypeListAt<2, Z0BC>;
  using HxZ0BC = TypeListAt<3, Z0BC>;
  using HyZ0BC = TypeListAt<4, Z0BC>;
  using HzZ0BC = TypeListAt<5, Z0BC>;

  using ExZ1BC = TypeListAt<0, Z1BC>;
  using EyZ1BC = TypeListAt<1, Z1BC>;
  using EzZ1BC = TypeListAt<2, Z1BC>;
  using HxZ1BC = TypeListAt<3, Z1BC>;
  using HyZ1BC = TypeListAt<4, Z1BC>;
  using HzZ1BC = TypeListAt<5, Z1BC>;

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
    DBG("updateE_bcs()");
    ExX0BC::apply();
    EyX0BC::apply();
    EzX0BC::apply();

    ExX1BC::apply();
    EyX1BC::apply();
    EzX1BC::apply();

    ExY0BC::apply();
    EyY0BC::apply();
    EzY0BC::apply();

    ExY1BC::apply();
    EyY1BC::apply();
    EzY1BC::apply();

    ExZ0BC::apply();
    EyZ0BC::apply();
    EzZ0BC::apply();

    ExZ1BC::apply();
    EyZ1BC::apply();
    EzZ1BC::apply();
  }

  static void updateH_bcs(auto& emdata, auto& bcdata) {
    DBG("updateH_bcs()");
    HxX0BC::apply();
    HyX0BC::apply();
    HzX0BC::apply();

    HxX1BC::apply();
    HyX1BC::apply();
    HzX1BC::apply();

    HxY0BC::apply();
    HyY0BC::apply();
    HzY0BC::apply();

    HxY1BC::apply();
    HyY1BC::apply();
    HzY1BC::apply();

    HxZ0BC::apply();
    HyZ0BC::apply();
    HzZ0BC::apply();

    HxZ1BC::apply();
    HyZ1BC::apply();
    HzZ1BC::apply();
  }


  static void advance(auto& emdata, auto& bcdata) {
    updateH(emdata);
    updateH_bcs(emdata, bcdata);

    updateE(emdata);
    updateE_bcs(emdata, bcdata);
  }
};

#endif //EM_SOLVER_H
