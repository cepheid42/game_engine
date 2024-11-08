//
// Created by cepheid on 10/17/24.
//

#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include "electromagnetics.param"
#include "em_traits.h"


//====== 1D Boundaries ========
//=============================
template<typename Array>
struct ReflectingBC {
  using array_t = Array;
  using value_t = typename Array::value_t;
  using dimension_t = typename Array::dimension_t;
  static constexpr size_t bc_depth = 0;
};

template<typename Array>
struct PeriodicBC {
  using array_t = Array;
  using value_t = typename Array::value_t;
  using dimension_t = typename Array::dimension_t;

  template<typename ARR>
  static void apply(ARR& f1, const size_t i) {
    DBG("PeriodicBC::apply()");
    const auto numInterior = f1.nx() - (2 * nHalo);
    const auto hi_idx = (f1.nx() - 1) - nHalo;
    const auto pm = i % numInterior;

    f1[nHalo - 1 - i] = f1[hi_idx - pm];
    f1[hi_idx + 1 + i] = f1[nHalo + pm];
  }

  template<typename ARR>
  requires is_empty_field<ARR, EmptyArray<typename ARR::value_t, ARR::dimension_t::value>>
  static void apply(ARR&, const size_t)
  { DBG("PeriodicBC::apply()::empty"); }
};

//template<typename Array, typename CurlA, typename CurlB, size_t DEPTH, typename... IDXS>
//struct PmlBC {
//  using array_t = Array;
//  using value_t = typename Array::value_t;
//  using dimension_t = typename Array::dimension_t;
//
//  static constexpr size_t bc_depth = DEPTH;
//
//  static void apply(auto& f1, size_t i) {
//    DBG("PMLBC::apply()");
//  }
//
//  template<typename ARR>
//  requires is_empty_field<ARR, EmptyArray<value_t, dimension_t::value>>
//  static void apply(ARR&, const size_t) // const auto&, const auto&, const auto&, const auto&, const auto&)
//  { DBG("PMLBC::apply()::empty"); }
//};


template<typename UpdateFunctor>
struct BCIntegrator1D {
  using update_t = UpdateFunctor;
  using array_t = typename update_t::array_t;
  using value_t = typename update_t::value_t;
  using dimension_t = typename update_t::dimension_t;

  // static constexpr size_t bc_depth = update_t::bc_depth;

  static void apply(auto& f1, const auto& o) {
    if constexpr (!std::same_as<UpdateFunctor, ReflectingBC<array_t>>) {
      DBG("BCIntegrator1D::apply()");
      for (size_t i = o.x0; i < o.x1; ++i) {
        update_t::apply(f1, i);
      }
    } else {
      DBG("BCIntegrator1D::apply()::empty");
    }
  }
};

//template<typename UpdateFunctor>
//struct BCIntegrator2D {
//  using update_t = UpdateFunctor;
//  using value_t = typename update_t::value_t;
//  using dimension_t = typename update_t::dimension_t;
//  using array_t = Array2D<value_t>;
//
//  static constexpr size_t bc_depth = UpdateFunctor::bc_depth;
//
//  static void apply(auto& f1) {//, const auto& d1, const auto& d2, const auto& c_d, auto& psi, const auto& b, const auto& c, const offset_t& o) {
//    if constexpr (!std::same_as<UpdateFunctor, ReflectingBC<array_t>>) {
//      DBG("BCIntegrator2D::apply()");
//      for (size_t i = f1.offsets.x0; i < f1.offsets.x1; ++i) {
//        for (size_t j = f1.offsets.y0; j < f1.offsets.y1; ++j) {
//          update_t::apply(f1, i, j);
//        }
//      }
//    } else {
//      DBG("BCIntegrator2D::apply()::empty");
//    }
//  }
//};
//
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
