//
// Created by cepheid on 10/17/24.
//

#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include "curl_operators.h"

template<typename T>
struct BCIntegratorNull {
  using value_t = typename T::value_t;
  using dimension_t = typename T::dimension_t;
  using array_t = EmptyArray<value_t, dimension_t::value>;

  static constexpr void apply(const auto&, const auto&, const auto&, const auto&, const auto&, const auto&, const auto&, const auto&) { DBG("BCIntegratorNull::apply()"); }
};

template<typename T, typename UpdateFunctor>
struct BCIntegrator1D {
  using value_t = typename T::value_t;
  using dimension_t = typename T::dimension_t;
  using array_t = Array1D<value_t>;
  using update_func = UpdateFunctor;

  static void apply(auto& f1, const auto& d1, const auto& d2, const auto& c_d, auto& psi, const auto& b, const auto& c, const auto& o) {
    DBG("BCIntegrator1D::apply()");
    for (size_t i = o.x0; i < o.x1; ++i) {
      update_func::apply(f1, d1, d2, c_d, psi, b, c, i);
    }
  }
};

template<typename T, typename UpdateFunctor>
struct BCIntegrator2D {
  using value_t = typename T::value_t;
  using dimension_t = typename T::dimension_t;
  using array_t = Array2D<value_t>;
  using update_func = UpdateFunctor;

  static void apply(auto& f1, const auto& d1, const auto& d2, const auto& c_d, auto& psi, const auto& b, const auto& c, const auto& o) {
    for (size_t i = o.x0; i < o.x1; ++i) {
      for (size_t j = o.y0; j < o.y1; ++j) {
        update_func::apply(f1, d1, d2, c_d, psi, b, c, i, j);
      }
    }
  }
};

//====== 1D Boundaries ========
//=============================
template<typename, bool=false, bool=false>
struct NoneBC {};

template<typename Array, typename... IDXS>
struct PeriodicBC {
  using value_t = typename Array::value_t;
  using dimension_t = typename Array::dimension_t;

  static void apply1D(auto& f1, const size_t p) {
    DBG("PeriodicBC::apply1D()");
    const auto numInterior = (f1.nx) - (2 * nHalo);
    const auto hi_idx = (f1.nx - 1) - (nHalo);

    // DBG(f1.nx, p, numInterior, nHalo, hi_idx);
    const auto pm = p % numInterior;

    // DBG(nHalo - 1 - p, hi_idx - pm, hi_idx + 1 + p, nHalo + pm);

    f1[nHalo - 1 - p] = f1[hi_idx - pm];
    f1[hi_idx + 1 + p] = f1[nHalo + pm];
  }

  static void apply2D(auto& f1, const size_t i, const size_t j) {
    // DBG("Periodic1D::apply1D()");
    // static auto numInterior = (f1.nx) - (2 * nHalo);
    // static auto lo_idx = nHalo;
    // static auto hi_idx = (f1.nx - 1) - (nHalo);
    //
    // for (size_t p = 0; p < nHalo; ++p) {
    //   // degenerate case
    //   const auto pm = p % numInterior;
    //   f1[lo_idx - 1 - p] = f1[hi_idx - pm];
    //   f1[hi_idx + 1 + p] = f1[lo_idx + pm];
    // }
  }


  static void apply3D(auto& f1, const size_t i, const size_t j, const size_t k) {}

  static void apply(auto& f1, const auto&, const auto&, const auto&, const auto&, const auto&, const auto&, IDXS... idxs) {
    DBG("PeriodicBC::apply()");
    if constexpr (dimension_t::value == 1) {
      apply1D(f1, idxs...);
    } else if constexpr (dimension_t::value == 2) {
      apply2D(f1, idxs...);
    } else {
      apply3D(f1, idxs...);
    }
  }

  static void apply(const auto&, const auto&, const auto&, const auto&, const auto&, const auto&, const auto&)
  requires is_empty_field<Array, EmptyArray<value_t, dimension_t::value>>
  { DBG("PeriodicBc::apply()::empty"); }
};

template<typename Array, bool HI, bool Forward, typename... IDXS>
struct PmlBC {
  // using array_t = Array;
  using value_t = typename Array::value_t;
  using dimension_t = typename Array::dimension_t;
  using curl1 = curl<Derivative::DX, Forward, size_t>;
  using curl2 = curl<Derivative::NoOp, Forward, size_t>;

  static constexpr bool hi_side = HI;
  static constexpr bool forward = Forward;

  static void apply1D(auto& f1, const auto& d1, const auto& d2, const auto& c_d, auto& psi, const auto& b, const auto& c, const size_t i) {
    // size_t start, stop;
    // if constexpr (HI) {
    //   start = (f1.nx - 1) - dPML;  // need -1?
    //   stop = f1.nx - 1;
    // } else {
    //   start = 0;
    //   stop = dPML;
    // }
    // const auto ipml = i - start;
    //
    // const auto self = b[ipml] * psi[ipml];
    // const auto diff1 = curl1::apply(d1, i);
    // const auto diff2 = curl2::apply(d2, i);
    // const auto diff = c[ipml] * (diff1 - diff2);
    //
    // psi[ipml] = self + diff;
    // f1[i] += c_d[i] * psi[ipml];
  }

  static void apply2D(auto& f1, const auto& d1, const auto& d2, const auto& c_d, auto& psi, const auto& b, const auto& c) {}
  static void apply3D(auto& f1, const auto& d1, const auto& d2, const auto& c_d, auto& psi, const auto& b, const auto& c) {}

  static void apply(auto& f1, const auto& d1, const auto& d2, const auto& c_d, auto& psi, const auto& b, const auto& c, IDXS... idxs) {
    if constexpr (dimension_t::value == 1) {
      apply1D(f1, d1, d2, c_d, psi, b, c, idxs...);
    } else if constexpr (dimension_t::value == 2) {
      apply2D(f1, d1, d2, c_d, psi, b, c, idxs...);
    } else {
      apply3D(f1, d1, d2, c_d, psi, b, c, idxs...);
    }
  }


  static void apply(auto&, const auto&, const auto&, const auto&, const auto&, const auto&)
  requires is_empty_field<Array, EmptyArray<value_t, dimension_t::value>>
  { DBG("PML1D::apply()::empty"); }
};



//====== 2D Boundaries ========
//=============================


//====== 3D Boundaries ========
//=============================

#endif //BOUNDARIES_H
