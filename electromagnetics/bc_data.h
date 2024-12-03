//
// Created by cepheid on 10/23/24.
//

#ifndef BC_DATA_H
#define BC_DATA_H

// #include "em_traits.h"
// #include "aydenstuff/array.h"
// #include <electromagnetics.h>
// #include "boundaries.h"

#include <array>

#include "electromagnetics.param"

using tf::types::Array1D;
using tf::types::Array2D;
using tf::types::Array3D;

template<typename T>
std::vector<T> linspace(T start, T stop, size_t n_points, const bool endpoint=true) {
  std::vector<T> result(n_points);
  if (endpoint) {
    n_points -= 1;
    result[result.size() - 1] = stop;
  }
  auto delta = (stop - start) / static_cast<T>(n_points);
  T val = start;
  for (size_t i = 0; i < n_points; ++i) {
    result[i] = val;
    val += delta;
  }
  return result;
}

template<typename Array>
struct NullData {
  using array_t = Array;
  using value_t = typename Array::value_t;
  using dimension_t = typename Array::dimension_t;

  explicit NullData(const auto&) {}
};


template<EMFace F>
size_t get_num_interior(const auto& f) {
  if constexpr (F == EMFace::X) {
    return f.nx() - (2 * nHalo);
  } else if constexpr (F == EMFace::Y) {
    return f.ny() - (2 * nHalo);
  } else {
    return f.nz() - (2 * nHalo);
  }
}

template<EMFace F>
size_t get_hi_index(const auto& f) {
  if constexpr (F == EMFace::X) {
    return f.nx() - 1 - nHalo;
  } else if constexpr (F == EMFace::Y) {
    return f.ny() - 1 - nHalo;
  } else {
    return f.nz() - 1 - nHalo;
  }
}

template<typename Array, EMFace F, EMSide S>
struct PeriodicData {
  using array_t = Array;
  using value_t = typename Array::value_t;
  using dimension_t = typename Array::dimension_t;

  explicit PeriodicData(const Array& f)
  : numInterior{get_num_interior<F>(f)},
    hi_idx{get_hi_index<F>(f)},
    offsets{get_offsets<F, S, nHalo>(f)}
  {}

  static constexpr size_t depth = nHalo;
  size_t numInterior;
  size_t hi_idx;
  IntegratorOffsets offsets;
};




template<typename Array, EMFace F, EMSide S, bool HF>
struct PMLData {
  using array_t = Array;
  using value_t = typename Array::value_t;
  using dimension_t = typename Array::dimension_t;

  explicit PMLData(const Array& f) requires (F == EMFace::X)
  : offsets{get_offsets<F, S, nPml>(f)},
    psi{nPml, f.ny(), f.nz()}
  {
    set_coefficients();
  }

  explicit PMLData(const Array& f) requires (F == EMFace::Y)
  : offsets{get_offsets<F, S, nPml>(f)},
    psi{f.nx(), nPml, f.nz()}
  {
    set_coefficients();
  }

  explicit PMLData(const Array& f) requires (F == EMFace::Z)
  : offsets{get_offsets<F, S, nPml>(f)},
    psi{f.nx(), f.ny(), nPml}
  {
    set_coefficients();
  }

  void set_coefficients() {
    auto d = linspace(1.0, 0.0, nPml, false);
    constexpr value_t hstep = 1.0 / (2.0 * static_cast<value_t>(nPml));

    if constexpr (HF) {
      for (auto& x: d) {
        x -= hstep;
      }
    }

    if constexpr (S == EMSide::Hi) {
      std::ranges::reverse(d);
    }

    constexpr auto c0 = 299792458.0;
    constexpr auto eps0 = 8.854187812813e-12;
    constexpr auto eta0 = 376.73031366686992;
    constexpr auto grade = 3.0;
    constexpr auto dx = 1.0 / 99.0;
    constexpr auto sigma_max = (0.8 * (grade + 1.0)) / (dx * eta0);
    constexpr auto alpha_max = 0.2;

    constexpr auto dt = cfl * dx / c0;

    // DBG(dt, dx, sigma_max, alpha_max);
    // DBG(d);

    std::vector<value_t> sigma_bc(d);
    std::vector<value_t> alpha_bc(d);
    std::vector<value_t> kappa_bc(d.size(), 1.0);

    for (auto& x: sigma_bc) {
      x = sigma_max * std::pow(x, grade);
    }

    for (auto& x: alpha_bc) {
      x = alpha_max * std::pow(1.0 - x, 1.0);
    }

    const auto coef1 = -dt / eps0;

    for (size_t i = 0; i < nPml; ++i) {
      b[i] = std::exp(coef1 * ((sigma_bc[i] / kappa_bc[i]) + alpha_bc[i]));
      c[i] = (sigma_bc[i] * (b[i] - 1.0)) / (dx * kappa_bc[i] * (sigma_bc[i] + (kappa_bc[i] * alpha_bc[i])));
    }

    // DBG(b, c);
  }

  static constexpr size_t depth = nPml;
  IntegratorOffsets offsets;
  array_t psi;
  std::array<value_t, nPml> b{};
  std::array<value_t, nPml> c{};
};


template<typename ex_t, typename ey_t, typename ez_t, typename hx_t, typename hy_t, typename hz_t>
struct FaceBCs {
  using value_t = typename ex_t::value_t;
  using dimension_t = typename ex_t::dimension_t;

  explicit FaceBCs(const auto& emdata)
  : Ex{emdata.Ex},
    Ey{emdata.Ey},
    Ez{emdata.Ez},
    Hx{emdata.Hx},
    Hy{emdata.Hy},
    Hz{emdata.Hz}
  {}

  ex_t Ex;
  ey_t Ey;
  ez_t Ez;
  hx_t Hx;
  hy_t Hy;
  hz_t Hz;
};

template<typename X0BC, typename X1BC, typename Y0BC, typename Y1BC, typename Z0BC, typename Z1BC>
struct BCData {
  using value_t = typename X0BC::value_t;
  using dimension_t = typename X0BC::dimension_t;

  explicit BCData(const auto& emdata)
  : x0{emdata},
    x1{emdata},
    y0{emdata},
    y1{emdata},
    z0{emdata},
    z1{emdata}
  {}

  X0BC x0;
  X1BC x1;
  Y0BC y0;
  Y1BC y1;
  Z0BC z0;
  Z1BC z1;
};

#endif //BC_DATA_H
