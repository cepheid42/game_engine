//
// Created by cepheid on 10/17/24.
//

#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include "curl_operators.h"
// #include "em_data.h"
// #include "aydenstuff/array.h"

constexpr size_t dPML = 10u;
constexpr size_t nHalo = 2u;

template<typename T, typename UpdateFunctor>
struct BCIntegrator1D {
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
struct BCIntegrator2D {
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

//====== 1D Boundaries ========
//=============================
template<typename Array, bool=false, bool=false>
struct NoneBC {
  using array_t = Array;
  using value_t = typename Array::value_t;
  using dimension_t = typename Array::dimension_t;

  static void apply(const auto&, const auto&, const auto&, const auto&, const auto&, const auto&) {}
};

template<typename Array, bool=false, bool=false>
struct PeriodicBC {
  // using array_t = Array;
  using value_t = typename Array::value_t;
  using dimension_t = typename Array::dimension_t;

  static void apply1D(auto& f1, const auto&, const auto&, const auto&, const auto&, const auto&) {
    DBG("Periodic1D::apply1D()");
    static auto numInterior = (f1.nx) - (2 * nHalo);
    static auto lo_idx = nHalo;
    static auto hi_idx = (f1.nx - 1) - (nHalo);

    for (size_t p = 0; p < nHalo; ++p) {
      // degenerate case
      // const auto pm = p % numInterior;
      // f1[lo_idx - 1 - p] = f1[hi_idx - pm];
      // f1[hi_idx + 1 + p] = f1[lo_idx + pm];
      f1[lo_idx - 1 - p] = f1[hi_idx - p];
      f1[hi_idx + 1 + p] = f1[lo_idx + p];
    }
  }
  static void apply2D(auto& f1, const auto& f2, auto& psi, const auto& c_f2, const auto& b, const auto& c) {}
  static void apply3D(auto& f1, const auto& f2, auto& psi, const auto& c_f2, const auto& b, const auto& c) {}

  static void apply(auto& f1, const auto& d1, const auto& d2, const auto& c_d, auto& psi, const auto& b, const auto& c) {
    DBG("Periodic1D::apply()");
    if constexpr (dimension_t::value == 1) {
      apply1D(f1);
    } else if constexpr (dimension_t::value == 2) {
      apply2D(f1);
    } else {
      apply3D(f1);
    }
  }

  static void apply(auto&, const auto&, const auto&, const auto&, const auto&, const auto&)
  requires is_empty_field<Array, EmptyArray<value_t, dimension_t::value>>
  { DBG("Periodic1D::apply()::empty"); }
};

template<typename Array, bool HI, bool Forward>
struct PmlBC {
  // using array_t = Array;
  using value_t = typename Array::value_t;
  using dimension_t = typename Array::dimension_t;
  using curl1 = curl<Derivative::DX, Forward, size_t>;
  using curl2 = curl<Derivative::NoOp, Forward, size_t>;

  static constexpr bool hi_side = HI;
  static constexpr bool forward = Forward;

  static void apply1D(auto& f1, const auto& d1, const auto& d2, const auto& c_d, auto& psi, const auto& b, const auto& c) {
    size_t start, stop;
    if constexpr (HI) {
      start = (f1.nx - 1) - dPML;  // need -1?
      stop = f1.nx - 1;
    } else {
      start = 0;
      stop = dPML;
    }

    for (size_t i = start; i < stop; ++i) {
      const auto ipml = i - start;

      const auto self = b[ipml] * psi[ipml];
      const auto diff1 = curl1::apply(d1, i);
      const auto diff2 = curl2::apply(d2, i);
      const auto diff = c[ipml] * (diff1 - diff2);

      psi[ipml] = self + diff;
      f1[i] += c_d[i] * psi[ipml];
    }
  }

  static void apply2D(auto& f1, const auto& d1, const auto& d2, const auto& c_d, auto& psi, const auto& b, const auto& c) {}
  static void apply3D(auto& f1, const auto& d1, const auto& d2, const auto& c_d, auto& psi, const auto& b, const auto& c) {}

  static void apply(auto& f1, const auto& d1, const auto& d2, const auto& c_d, auto& psi, const auto& b, const auto& c) {
    if constexpr (dimension_t::value == 1) {
      apply1D(f1, d1, d2, c_d, psi, b, c);
    } else if constexpr (dimension_t::value == 2) {
      apply2D(f1, d1, d2, c_d, psi, b, c);
    } else {
      apply3D(f1, d1, d2, c_d, psi, b, c);
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
