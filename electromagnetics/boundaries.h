//
// Created by cepheid on 10/17/24.
//

#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include <string>

#include "electromagnetics.param"
// #include "em_traits.h"


//====== 1D Boundaries ========
//=============================
// template<typename Array>
// struct ReflectingBC {
//   using array_t = Array;
//   using value_t = typename Array::value_t;
//   using dimension_t = typename Array::dimension_t;
// };

template<typename... IDXS>
struct PeriodicBC {

  static void apply(const std::string& msg) {
    DBG("PeriodicBC::apply(" + msg + ")");
    // const auto numInterior = f1.nx() - (2 * nHalo);
    // const auto hi_idx = (f1.nx() - 1) - nHalo;
    // const auto pm = i % numInterior;
    //
    // f1[nHalo - 1 - i] = f1[hi_idx - pm];
    // f1[hi_idx + 1 + i] = f1[nHalo + pm];
  }

};

template<Derivative CURL1, Derivative CURL2, bool Forward, typename... IDXS>
struct PmlBC {
  using CurlA = curl<CURL1, Forward, IDXS...>;
  using CurlB = curl<CURL2, Forward, IDXS...>;

  static void apply(auto& bc, auto& f1, auto& f2, auto& c1, size_t i, size_t j, size_t)
  requires (sizeof...(IDXS) == 2)
  {
    DBG("PMLBC_2D::apply()");
    const auto ipml = i - start;

    const auto self = bc.b[ipml] * bc.psi(ipml, j);
    const auto diff1 = CurlA::apply(f2, i );
    const auto diff2 = CurlB::apply(f2, i);
    const auto diff = bc.c[ipml] * (diff1 - diff2);
  }

};


template<typename T>
struct BCIntegratorNull {
  using value_t = typename T::value_t;
  using dimension_t = typename T::dimension_t;
  using array_t = EmptyArray<value_t, dimension_t::value>;

  static constexpr void apply() { DBG("BCIntegratorNull::apply()"); }
};

// template<typename T, typename UpdateFunctor>
// struct BCIntegrator1D {
//   using update_t = UpdateFunctor;
//   using value_t = typename T::value_t;
//   using dimension_t = typename T::dimension_t;
//   using array_t = Array1D<value_t>;
//
//   static void apply() {
//     DBG("BCIntegrator1D::apply()");
//     // for (size_t i = o.x0; i < o.x1; ++i) {
//     //   update_t::apply(f1, i);
//     // }
//   }
// };
//
// template<typename T, typename UpdateFunctor>
// struct BCIntegrator2D {
//   using update_t = UpdateFunctor;
//   using value_t = typename T::value_t;
//   using dimension_t = typename T::dimension_t;
//   using array_t = Array2D<value_t>;
//
//   static void apply() {
//     DBG("BCIntegrator2D::apply()" );
//     // for (size_t i = f1.offsets.x0; i < f1.offsets.x1; ++i) {
//     //   for (size_t j = f1.offsets.y0; j < f1.offsets.y1; ++j) {
//     //     update_t::apply(f1, i, j);
//     //   }
//     // }
//   }
// };

template<typename T, typename UpdateFunctor>
struct BCIntegrator {
  using update_t = UpdateFunctor;
  using value_t = typename T::value_t;
  using dimension_t = typename T::dimension_t;
  using array_t = typename T::array_t;

  static void apply(auto& bc, auto& f1, auto& f2, auto& c1) {
    DBG("BCIntegrator3D::apply()");
    for (size_t i = bc.offsets.x0; i < bc.offsets.x1; ++i) {
      for (size_t j = bc.offsets.y0; j < bc.offsets.y1; ++j) {
        for (size_t k = bc.offsets.z0; k < bc.offsets.z1; ++k) {
          update_t::apply(bc, f1, f2, c1, i, j, k);
        }
      }
    }
  }
};


template<EMFace F, typename ex_t, typename ey_t, typename ez_t, typename hx_t, typename hy_t, typename hz_t>
struct EMBoundary {
  using value_t = typename ex_t::value_t;
  using dimension_t = typename ex_t::dimension_t;
  using empty_t = EmptyArray<value_t, dimension_t::value>;

  static constexpr empty_t empty{};

  static void updateE(auto& emdata, auto& bcdata)
  requires (F == EMFace::X)
  {
    ex_t::apply(bcdata.Ex, emdata.Ex, empty, empty);
    ey_t::apply(bcdata.Ey, emdata.Ey, emdata.Hz, emdata.Ceyh);
    ez_t::apply(bcdata.Ez, emdata.Ez, emdata.Hy, emdata.Cezh);
  }

  static void updateE(auto& emdata, auto& bcdata)
  requires (F == EMFace::Y)
  {
    ex_t::apply(bcdata.Ex, emdata.Ex, emdata.Hz, emdata.Cexh);
    ey_t::apply(bcdata.Ey, emdata.Ey, empty, empty);
    ez_t::apply(bcdata.Ez, emdata.Ez, emdata.Hx, emdata.Cezh);
  }

  static void updateE(auto& emdata, auto& bcdata)
  requires (F == EMFace::Z)
  {
    ex_t::apply(bcdata.Ex, emdata.Ex, emdata.Hz, emdata.Cexh);
    ey_t::apply(bcdata.Ey, emdata.Ey, emdata.Hz, emdata.Ceyh);
    ez_t::apply(bcdata.Ez, emdata.Ez, empty, empty);
  }

  static void updateH(auto& emdata, auto& bcdata) {}
};

#endif //BOUNDARIES_H
