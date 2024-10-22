//
// Created by cepheid on 10/17/24.
//

#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include <electromagnetics.h>
#include <electromagnetics.h>

#include "curl_operators.h"
#include "em_data.h"
#include "aydenstuff/array.h"

constexpr size_t dPML = 10u;
constexpr size_t nHalo = 2u;

//====== 1D Boundaries ========
//=============================
template<typename Array, bool=false, bool=false>
struct None1D {
  using array_t = Array;
  using value_t = typename Array::value_t;
  using dimension_t = typename Array::dimension_t;

  static void apply(const auto&, const auto&, const auto&, const auto&, const auto&, const auto&) { DBG("None1D::apply()"); }
};

template<typename Array, bool=false, bool=false>
struct Periodic1D {
  using array_t = Array;
  using value_t = typename Array::value_t;
  using dimension_t = typename Array::dimension_t;


  static void apply(array_t& f1, const auto&, const auto&, const auto&, const auto&, const auto&) {
    DBG("Periodic1D::apply()");
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

  static void apply(array_t&, const auto&, const auto&, const auto&, const auto&, const auto&)
  requires is_empty_field<array_t, EmptyArray<value_t, dimension_t::value>>
  { DBG("Periodic1D::apply()::empty"); }
};

template<typename Array, bool Forward, bool Subtract>
struct PML1D {
  using array_t = Array;
  using value_t = typename Array::value_t;
  using dimension_t = typename Array::dimension_t;

  using curl1 = curl<Derivative::DX, Forward, size_t>;

  static void apply(array_t& f1, const array_t& f2, array_t& psi, const array_t& c_f2, const array_t& b, const array_t& c) {
    for (size_t i = 0; i < dPML; ++i) {
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

  static void apply(array_t&, const auto&, const auto&, const auto&, const auto&, const auto&)
  requires is_empty_field<array_t, EmptyArray<value_t, dimension_t::value>>
  { DBG("PML1D::apply()::empty"); }
};

enum class BoundarySide { X, Y, Z };

template<typename BCType, BoundarySide S, FieldComponent EXF, FieldComponent EYF, FieldComponent EZF, FieldComponent HXF, FieldComponent HYF, FieldComponent HZF>
struct Boundary {
  using value_t = typename EXF::array_t::value_t;
  using dimension_t = typename EXF::array_t::dimension_t;
  using empty_t = EmptyArray<value_t, dimension_t::value>;

  static constexpr empty_t empty{};

  using ex_t = typename EXF::array_t;
  using ey_t = typename EYF::array_t;
  using ez_t = typename EZF::array_t;
  using hx_t = typename HXF::array_t;
  using hy_t = typename HYF::array_t;
  using hz_t = typename HZF::array_t;

  template<BoundarySide SIDE>
  explicit Boundary(auto& emdata)
  : ex(&emdata.Ex), cexh(&emdata.Cexh),
    ey(&emdata.Ey), ceyh(&emdata.Ceyh),
    ez(&emdata.Ez), cezh(&emdata.Cezh),
    hx(&emdata.Hx), chxe(&emdata.Chxe),
    hy(&emdata.Hy), chye(&emdata.Chye),
    hz(&emdata.Hz), chze(&emdata.Chze)
  {
    if constexpr (!std::is_same_v<BCType, Periodic1D<Array1D<value_t>>> and !std::is_same_v<BCType, None1D<Array1D<value_t>>>) {
      if constexpr (!std::is_same_v<ex_t, EmptyArray<value_t, dimension_t::value>>) {
        // psi_ex = ex_t{ex.nx};
        // bex = ex_t{ex.nx};
        // cex = ex_t{ex.nx};
        // todo: this is not correct for the X-face stuff. Also, its getting way too complicated.
      }
    }

  }

  void updateNone() {}
  
  template<BoundarySide SIDE>
  void updateE_periodic() {}

  template<BoundarySide SIDE>
  void updateE_pml() {
    if constexpr (SIDE == BoundarySide::X) {
      BCType::apply(&ey, &hz, &psi_ey, &ceyh, &bey, &cey);
      BCType::apply(&ez, &hy, &psi_ez, &cezh, &bez, &cez);
    } else if constexpr (SIDE == BoundarySide::Y) {
      BCType::apply(&ex, &hz, &psi_ex, &cexh, &bex, &cex);
      BCType::apply(&ez, &hx, &psi_ez, &cezh, &bez, &cez);
    } else {
      BCType::apply(&ex, &hy, &psi_ex, &cexh, &bex, &cex);
      BCType::apply(&ey, &hx, &psi_ey, &ceyh, &bez, &cez);
    }
  }
  
  void updateE() {
    if constexpr (std::is_same_v<BCType, Periodic1D<Array1D<value_t>>>) {
      updateE_periodic<S>();
    } else if constexpr (std::is_same_v<BCType, None1D<Array1D<value_t>>>) {
      updateNone();
    } else {
      updateE_pml<S>();
    }
  }
  void updateH() {}

  ex_t* ex;
  ex_t* cexh;
  
  ey_t* ey;
  ey_t* ceyh;
  
  ez_t* ez;
  ez_t* cezh;
  
  hx_t* hx;
  hx_t* chxe;
  
  hy_t* hy;
  hy_t* chye;
  
  hz_t* hz;
  hz_t* chze;
  
  ex_t psi_ex, bex, cex;
  ey_t psi_ey, bey, cey;
  ez_t psi_ez, bez, cez;

  hx_t psi_hx, bhx, chx;
  hy_t psi_hy, bhy, chy;
  hz_t psi_hz, bhz, chz;
};

//====== 2D Boundaries ========
//=============================


//====== 3D Boundaries ========
//=============================

#endif //BOUNDARIES_H
