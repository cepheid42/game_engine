//
// Created by cepheid on 10/17/24.
//

#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include "curl_operators.h"
#include "em_data.h"
#include "aydenstuff/array.h"

constexpr size_t dPML = 10;


//====== 1D Boundaries ========
//=============================
template<typename Array, bool, bool>
struct Periodic1D {
  using array_t = Array;
  using value_t = typename Array::value_t;

  static void apply(array_t& f1, const auto&, const auto&, const auto&, const auto&, const auto&, size_t nHalo) {
    static auto numInterior = (f1.nx) - (2 * nHalo);
    static auto lo_idx = nHalo;
    static auto hi_idx = (f1.nx - 1) - (nHalo);

    // DBG(numInterior, lo_idx, hi_idx);
    for (size_t p = 0; p < nHalo; ++p) {
      // const auto pm = p % numInterior;
      // DBG(p);
      // DBG(lo_idx - 1 - p, hi_idx - p);
      // DBG(hi_idx + 1 + p, lo_idx + p);

      f1[lo_idx - 1 - p] = f1[hi_idx - p];
      f1[hi_idx + 1 + p] = f1[lo_idx + p];
    }
  }
};

template<typename Array, bool Forward, bool Subtract>
struct PML1D {
  using array_t = Array;
  using value_t = typename Array::value_t;

  using curl1 = curl<Derivative::DX, Forward, size_t>;

  static void apply(array_t& f1, array_t& f2, array_t& psi, const array_t& c_f2, const array_t& b, const array_t& c, size_t depth) {
    for (size_t i = 0; i < depth; ++i) {
      const auto self = b[i] * psi[i];
      const auto curl = c[i] * curl1::apply(f2, i);
      psi[i] = self + curl;

      if constexpr (Subtract) {
        f1[i] -= c_f2[i] * psi[i];
      } else {
        f1[i] += c_f2[i] * psi[i];
      }
    }
  }
};




//====== 2D Boundaries ========
//=============================


//====== 3D Boundaries ========
//=============================

#endif //BOUNDARIES_H
