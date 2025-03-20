#ifndef EM_BOUNDARIES_HPP
#define EM_BOUNDARIES_HPP

#include "program_params.hpp"
#include "em_params.hpp"
#include "array.hpp"
#include "constants.hpp"
#include "math_utils.hpp"

#include <vector>
#include <algorithm>
#include <array>

namespace tf::electromagnetics {
  enum class EMFace { X, Y, Z};
  enum class EMSide { Lo, Hi };

  template<EMSide S>
  std::array<std::size_t, 6> get_x_offsets(const std::size_t nx, const std::size_t ny, const std::size_t nz) {
    if constexpr (S == EMSide::Lo) {
      return {0, PMLDepth, 0, ny, 0, nz};
    } else {
      return {nx - PMLDepth - 1, nx - 1, 0, ny, 0, nz};
    }
  }

  template<EMSide S>
  std::array<std::size_t, 6> get_y_offsets(const std::size_t nx, const std::size_t ny, const std::size_t nz) {
    if constexpr (S == EMSide::Lo) {
      return {0, nx, 0, PMLDepth, 0, nz};
    } else {
      return {0, nx, ny - PMLDepth - 1, ny - 1, 0, nz};
    }
  }

  template<EMSide S>
  std::array<std::size_t, 6> get_z_offsets(const std::size_t nx, const std::size_t ny, const std::size_t nz) {
    if constexpr (S == EMSide::Lo) {
      return {0, nx, 0, ny, 0, PMLDepth};
    } else {
      return {0, nx, 0, ny, nz - PMLDepth - 1, nz - 1};
    }
  }

  template<EMFace F, EMSide S, bool isE>
  auto getOffsets(const auto& f) {
    std::array<std::size_t, 6> result{};

    if constexpr (F == EMFace::X) {
      result = get_x_offsets<S>(f.nx(), f.ny(), f.nz());
    } else if constexpr (F == EMFace::Y) {
      result = get_y_offsets<S>(f.nx(), f.ny(), f.nz());
    } else {
      result = get_z_offsets<S>(f.nx(), f.ny(), f.nz());
    }

    constexpr auto lo = static_cast<std::size_t>(isE == (S == EMSide::Lo));
    constexpr auto hi = static_cast<std::size_t>(isE != (S == EMSide::Lo));

    std::size_t lo_idx, hi_idx;
    if constexpr (F == EMFace::X) {
      lo_idx = 0; // x0
      hi_idx = 1; // x1
    } else if constexpr (F == EMFace::Y) {
      lo_idx = 2; // y0
      hi_idx = 3; // y1
    } else {
      lo_idx = 4; // z0
      hi_idx = 5; // z1
    }

    result[lo_idx] += lo;
    result[hi_idx] -= hi;

    return result;
  }


  template<EMFace F, EMSide S>
  struct PMLData {
    using offset_t = std::array<std::size_t, 6>;
    using coeffs_t = std::array<double, PMLDepth>;

    explicit PMLData() = delete;

    PMLData(const std::size_t nx_, const std::size_t ny_, const std::size_t nz_, const offset_t& offsets_, const bool isE)
    : psi(nx_, ny_, nz_),
      offsets(offsets_)
    {
      init_coefficients(isE);
    }

    PMLData(const Array3D<double>& f, const offset_t& offset, const bool isE)
    requires (F == EMFace::X)
    : PMLData(PMLDepth, f.ny(), f.nz(), offset, isE)
    {}

    PMLData(const Array3D<double>& f, const offset_t& offset, const bool isE)
    requires (F == EMFace::Y)
    : PMLData(f.nx(), PMLDepth, f.nz(), offset, isE)
    {}

    PMLData(const Array3D<double>& f, const offset_t& offset, const bool isE)
    requires (F == EMFace::Z)
    : PMLData(f.nx(), f.ny(), PMLDepth, offset, isE)
    {}

    void init_coefficients(const bool isE) {
      std::vector<double> d = math::linspace(1.0, 0.0, PMLDepth, false);

      if (!isE) {
        constexpr auto hstep = 1.0 / (2.0 * static_cast<double>(PMLDepth));
        for (auto& x: d) { x -= hstep; }
      }

      constexpr auto sigma_max = (0.8 * (PMLGrade + 1.0)) / (dx * constants::eta0);

      const auto sigma_d = calculate_sigma(d, sigma_max);
      const auto alpha_d = calculate_alpha(d);
      calculate_coeffs(sigma_d, alpha_d);

      if constexpr (S == EMSide::Hi) {
        std::ranges::reverse(b);
        std::ranges::reverse(c);
      }
    }

    static std::vector<double> calculate_sigma(const std::vector<double>& d, const double sigma_max) {
      auto sigma_bc(d);
      for (auto& x: sigma_bc) {
        x = sigma_max * std::pow(x, PMLGrade);
      }
      return sigma_bc;
    }

    static std::vector<double> calculate_alpha(const std::vector<double>& d) {
      auto alpha_bc(d);
      for (auto& x: alpha_bc) {
        x = PMLAlphaMax * std::pow(1.0 - x, 1.0);
      }
      return alpha_bc;
    }

    void calculate_coeffs(const std::vector<double>& sigma, const std::vector<double>& alpha) {
      constexpr auto coef1 = -dt / constants::eps0;

      for (auto i = 0zu; i < PMLDepth; i++) {
        constexpr auto kappa_bc = 1.0;
        b[i] = std::exp(coef1 * ((sigma[i] / kappa_bc) + alpha[i]));
        c[i] = (sigma[i] * (b[i] - 1.0)) / (kappa_bc * (sigma[i] + (kappa_bc * alpha[i])));
      }
    }

    Array3D<double> psi;
    offset_t offsets;
    coeffs_t b{};
    coeffs_t c{};
  };

  template<EMFace F, EMSide S>
  struct FaceBC {
    // template<bool isE> using getOffsets = get_pml_offsets<F, isE>;

    explicit FaceBC() = delete;

    explicit FaceBC(const auto& emdata)
    : Ex{emdata.Ex, getOffsets<F, S, true>(emdata.Ex), true},
      Ey{emdata.Ey, getOffsets<F, S, true>(emdata.Ey), true},
      Ez{emdata.Ez, getOffsets<F, S, true>(emdata.Ez), true},
      Hx{emdata.Hx, getOffsets<F, S, false>(emdata.Hx), false},
      Hy{emdata.Hy, getOffsets<F, S, false>(emdata.Hy), false},
      Hz{emdata.Hz, getOffsets<F, S, false>(emdata.Hz), false}
    {}

    PMLData<F, S> Ex;
    PMLData<F, S> Ey;
    PMLData<F, S> Ez;
    PMLData<F, S> Hx;
    PMLData<F, S> Hy;
    PMLData<F, S> Hz;
  };


  struct BCData {
    explicit BCData(const auto& emdata)
    : x0(emdata),
      y0(emdata),
      z0(emdata),
      x1(emdata),
      y1(emdata),
      z1(emdata)
    {}

    FaceBC<EMFace::X, EMSide::Lo> x0;
    FaceBC<EMFace::Y, EMSide::Lo> y0;
    FaceBC<EMFace::Z, EMSide::Lo> z0;
    FaceBC<EMFace::X, EMSide::Hi> x1;
    FaceBC<EMFace::Y, EMSide::Hi> y1;
    FaceBC<EMFace::Z, EMSide::Hi> z1;
  };
}

#endif //EM_BOUNDARIES_HPP
