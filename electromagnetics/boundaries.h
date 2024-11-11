//
// Created by cepheid on 10/17/24.
//

#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include "electromagnetics.param"
#include "em_traits.h"


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

  template<typename ARR>
  static void apply(ARR& f1, const size_t i) {
    DBG("PeriodicBC::apply()");
    const auto numInterior = f1.nx() - (2 * nHalo);
    const auto hi_idx = (f1.nx() - 1) - nHalo;
    const auto pm = i % numInterior;

    f1[nHalo - 1 - i] = f1[hi_idx - pm];
    f1[hi_idx + 1 + i] = f1[nHalo + pm];
  }

};

template<Derivative CURL1, Derivative CURL2, bool Forward, typename... IDXS>
struct PmlBC {
  using CurlA = curl<CURL1, Forward, IDXS...>;
  using CurlB = curl<CURL2, Forward, IDXS...>;

  static void apply(auto& f1, size_t i) {
    DBG("PMLBC::apply()");
  }

};


template<typename T>
struct BCIntegratorNull {
  using value_t = typename T::value_t;
  using dimension_t = typename T::dimension_t;
  using array_t = EmptyArray<value_t, dimension_t::value>;

  static constexpr void apply() {}
};

template<typename T, typename UpdateFunctor>
struct BCIntegrator1D {
  using update_t = UpdateFunctor;
  using value_t = typename T::value_t;
  using dimension_t = typename T::dimension_t;
  using array_t = Array1D<value_t>;

  // static constexpr size_t bc_depth = update_t::bc_depth;

  static void apply(auto& f1, const auto& o) {
    DBG("BCIntegrator1D::apply()");
    for (size_t i = o.x0; i < o.x1; ++i) {
      update_t::apply(f1, i);
    }
  }
};

template<typename UpdateFunctor>
struct BCIntegrator2D {
  using update_t = UpdateFunctor;
  using value_t = typename update_t::value_t;
  using dimension_t = typename update_t::dimension_t;
  using array_t = Array2D<value_t>;

  static constexpr size_t bc_depth = UpdateFunctor::bc_depth;

  static void apply(auto& f1) {//, const auto& d1, const auto& d2, const auto& c_d, auto& psi, const auto& b, const auto& c, const offset_t& o) {
    DBG("BCIntegrator2D::apply()");
    for (size_t i = f1.offsets.x0; i < f1.offsets.x1; ++i) {
      for (size_t j = f1.offsets.y0; j < f1.offsets.y1; ++j) {
        update_t::apply(f1, i, j);
      }
    }
  }
};

//template<typename UpdateFunctor>
//struct BCIntegrator3D {
//  using update_t = UpdateFunctor;
//  using value_t = typename update_t::value_t;
//  using dimension_t = typename update_t::dimension_t;
//  using array_t = Array3D<value_t>;
//
//  static constexpr size_t bc_depth = UpdateFunctor::bc_depth;
//
//  static void apply(auto& f1) { //, const auto& d1, const auto& d2, const auto& c_d, auto& psi, const auto& b, const auto& c, const offset_t& o) {
//    if constexpr (!std::same_as<UpdateFunctor, ReflectingBC<array_t>>) {
//      DBG("BCIntegrator3D::apply()");
//      for (size_t i = f1.offsets.x0; i < f1.offsets.x1; ++i) {
//        for (size_t j = f1.offsets.y0; j < f1.offsets.y1; ++j) {
//          for (size_t k = f1.offsets.z0; k < f1.offsets.z1; ++k) {
//            update_t::apply(f1, i, j, k);
//          }
//        }
//      }
//    } else {
//      DBG("BCIntegrator3D::apply()::empty");
//    }
//  }
//};

#endif //BOUNDARIES_H
