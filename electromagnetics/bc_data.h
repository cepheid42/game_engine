//
// Created by cepheid on 10/23/24.
//

#ifndef BC_DATA_H
#define BC_DATA_H

#include <array>

#include "electromagnetics.param"

namespace tf::electromagnetics::boundaries
{
  // todo: This will exist in the utility.h file
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

    explicit NullData(const auto&...) {}
  };

  template<typename Array, EMFace F, EMSide S>
  struct PeriodicData {
    using array_t = Array;
    using value_t = typename Array::value_t;
    using dimension_t = typename Array::dimension_t;
    using offset_t = typename tf::electromagnetics::types::IntegratorOffsets;

    explicit PeriodicData(const Array& f, const value_t, const value_t)
    : numInterior{get_num_interior<F>(f)},
      hi_idx{get_hi_index<F>(f)},
      offsets{get_offsets<F, S, nHalo>(f)}
    {}

    static size_t get_num_interior(const auto& f) {
      if constexpr (F == EMFace::X) {
        return f.nx() - (2 * nHalo);
      } else if constexpr (F == EMFace::Y) {
        return f.ny() - (2 * nHalo);
      } else {
        return f.nz() - (2 * nHalo);
      }
    }

    static size_t get_hi_index(const auto& f) {
      if constexpr (F == EMFace::X) {
        return f.nx() - 1 - nHalo;
      } else if constexpr (F == EMFace::Y) {
        return f.ny() - 1 - nHalo;
      } else {
        return f.nz() - 1 - nHalo;
      }
    }

    size_t numInterior;
    size_t hi_idx;
    offset_t offsets;
  };


  template<typename Array, EMFace F, EMSide S, bool HField>
  struct PMLData {
    using array_t = Array;
    using value_t = typename Array::value_t;
    using dimension_t = typename Array::dimension_t;
    using offset_t = typename tf::electromagnetics::types::IntegratorOffsets;

    explicit PMLData(const Array& f, const value_t dt, const value_t dx) requires (F == EMFace::X)
    : offsets{tf::electromagnetics::types::get_offsets<F, S, nPml>(f)},
      psi{nPml, f.ny(), f.nz()}
    {
      set_coefficients(dt, dx);
    }

    explicit PMLData(const Array& f, const value_t dt, const value_t dx) requires (F == EMFace::Y)
    : offsets{tf::electromagnetics::types::get_offsets<F, S, nPml>(f)},
      psi{f.nx(), nPml, f.nz()}
    {
      set_coefficients(dt, dx);
    }

    explicit PMLData(const Array& f, const value_t dt, const value_t dx) requires (F == EMFace::Z)
    : offsets{tf::electromagnetics::types::get_offsets<F, S, nPml>(f)},
      psi{f.nx(), f.ny(), nPml}
    {
      set_coefficients(dt, dx);
    }

    void set_coefficients(const value_t dt, const value_t dx) {
      // todo: These will all be in the constants header
      constexpr auto eps0 = 8.854187812813e-12;
      constexpr auto eta0 = 376.73031366686992;

      auto d = linspace(1.0, 0.0, nPml, false);
      auto coef1 = -dt / eps0;

      if constexpr (HField) {
        coef1 /= 2.0; // H field update is split into two steps, to need half dt
        constexpr value_t hstep = 1.0 / (2.0 * static_cast<value_t>(nPml));
        for (auto& x: d) {
          x -= hstep;
        }
      }

      if constexpr (S == EMSide::Hi) {
        std::ranges::reverse(d);
      }


      // todo: Add relative Mu and Eps to this to make it work for materials
      //       e.g. (0.8 * (grade + 1.0)) / (dx[i] * eta0 * sqrt(eps_r[i] * mu_r[i]))
      const auto sigma_max = (0.8 * (PMLGrade + 1.0)) / (dx * eta0);

      std::vector<value_t> sigma_bc(d);
      std::vector<value_t> alpha_bc(d);
      std::vector<value_t> kappa_bc(d.size(), 1.0);

      for (auto& x: sigma_bc) {
        x = sigma_max * std::pow(x, PMLGrade);
      }

      // todo: Kappa is hard coded to 1 for now, since it's not clear how useful it is to do otherwise
      //       and it couples the Solver and BC's together in an annoying way.
      // for (auto& x: kappa_bc) {
      //   x = 1.0 + (PMLKappaMax - 1.0) * std::pow(x, PMLGrade);
      // }

      for (auto& x: alpha_bc) {
        x = PMLAlphaMax * std::pow(1.0 - x, 1.0);
      }

      // todo: Pulling the dx term out of the 'c' coefficients to add it into the global coefficients (Cexhy, etc...)
      //       should require adding the dx term to the 'b' coefficients. However, this breaks the PML's entirely,
      //       so don't do that.
      //       e.g. This works and I don't know why. Do not touch.
      for (size_t i = 0; i < nPml; ++i) {
        b[i] = std::exp(coef1 * ((sigma_bc[i] / kappa_bc[i]) + alpha_bc[i]));
        c[i] = (sigma_bc[i] * (b[i] - 1.0)) / (kappa_bc[i] * (sigma_bc[i] + (kappa_bc[i] * alpha_bc[i])));
      }
    }

    offset_t offsets;
    array_t psi;
    std::array<value_t, nPml> b{};
    std::array<value_t, nPml> c{};
  };


  template<typename ex_t, typename ey_t, typename ez_t, typename hx_t, typename hy_t, typename hz_t>
  struct FaceBCs {
    using value_t = typename ex_t::value_t;
    using dimension_t = typename ex_t::dimension_t;

    // todo: dt/dx will probably get replaced with mesh object or something.
    explicit FaceBCs(const auto& emdata, const value_t dt, const value_t dx)
    : Ex{emdata.Ex, dt, dx},
      Ey{emdata.Ey, dt, dx},
      Ez{emdata.Ez, dt, dx},
      Hx{emdata.Hx, dt, dx},
      Hy{emdata.Hy, dt, dx},
      Hz{emdata.Hz, dt, dx}
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

    explicit BCData(const auto& emdata, const value_t dt, const value_t dx)
    : x0{emdata, dt, dx},
      x1{emdata, dt, dx},
      y0{emdata, dt, dx},
      y1{emdata, dt, dx},
      z0{emdata, dt, dx},
      z1{emdata, dt, dx}
    {}

    X0BC x0;
    X1BC x1;
    Y0BC y0;
    Y1BC y1;
    Z0BC z0;
    Z1BC z1;
  };
}
#endif //BC_DATA_H
