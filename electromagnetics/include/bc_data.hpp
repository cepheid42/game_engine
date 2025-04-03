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

  template<EMSide S, std::size_t Depth>
  std::array<std::size_t, 6> get_x_offsets(const std::size_t nx, const std::size_t ny, const std::size_t nz) {
    if constexpr (S == EMSide::Lo) {
      return {0, Depth, 0, ny, 0, nz};
    } else {
      return {nx - Depth - 1, nx - 1, 0, ny, 0, nz};
    }
  }

  template<EMSide S, std::size_t Depth>
  std::array<std::size_t, 6> get_y_offsets(const std::size_t nx, const std::size_t ny, const std::size_t nz) {
    if constexpr (S == EMSide::Lo) {
      return {0, nx, 0, Depth, 0, nz};
    } else {
      return {0, nx, ny - Depth - 1, ny - 1, 0, nz};
    }
  }

  template<EMSide S, std::size_t Depth>
  std::array<std::size_t, 6> get_z_offsets(const std::size_t nx, const std::size_t ny, const std::size_t nz) {
    if constexpr (S == EMSide::Lo) {
      return {0, nx, 0, ny, 0, Depth};
    } else {
      return {0, nx, 0, ny, nz - Depth - 1, nz - 1};
    }
  }

  template<EMFace F, EMSide S, std::size_t Depth>
  auto getOffsets(const auto& f) {
    if constexpr (F == EMFace::X) {
      return get_x_offsets<S, Depth>(f.nx(), f.ny(), f.nz());
    } else if constexpr (F == EMFace::Y) {
      return get_y_offsets<S, Depth>(f.nx(), f.ny(), f.nz());
    } else {
      return get_z_offsets<S, Depth>(f.nx(), f.ny(), f.nz());
    }
  }

  template<EMFace F, EMSide S, bool isE>
  auto getPMLOffsets(const auto& f) {
    auto result = getOffsets<F, S, PMLDepth>(f);
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
    using coeffs_t = std::array<compute_t, PMLDepth>;

    explicit PMLData() = delete;

    PMLData(const std::size_t nx_, const std::size_t ny_, const std::size_t nz_, const offset_t& offsets_, const bool isE)
    : psi(nx_, ny_, nz_),
      offsets(offsets_)
    {
      init_coefficients(isE);
    }

    PMLData(const Array3D<compute_t>& f, const offset_t& offset, const bool isE)
    requires (F == EMFace::X)
    : PMLData(PMLDepth, f.ny(), f.nz(), offset, isE)
    {}

    PMLData(const Array3D<compute_t>& f, const offset_t& offset, const bool isE)
    requires (F == EMFace::Y)
    : PMLData(f.nx(), PMLDepth, f.nz(), offset, isE)
    {}

    PMLData(const Array3D<compute_t>& f, const offset_t& offset, const bool isE)
    requires (F == EMFace::Z)
    : PMLData(f.nx(), f.ny(), PMLDepth, offset, isE)
    {}

    void init_coefficients(const bool isE) {
      std::vector<compute_t> d = math::linspace(1.0_fp, 0.0_fp, PMLDepth, false);

      if (!isE) {
        constexpr auto hstep = 1.0_fp / (2.0_fp * static_cast<compute_t>(PMLDepth));
        for (auto& x: d) { x -= hstep; }
      }

      constexpr auto eta0 = static_cast<compute_t>(constants::eta0);
      constexpr auto pml_grade = static_cast<compute_t>(PMLGrade);
      constexpr auto sigma_max = (0.8_fp * (pml_grade + 1.0_fp)) / (dx * eta0);

      const auto sigma_d = calculate_sigma(d, sigma_max);
      const auto alpha_d = calculate_alpha(d);
      calculate_coeffs(sigma_d, alpha_d);

      if constexpr (S == EMSide::Hi) {
        std::ranges::reverse(b);
        std::ranges::reverse(c);
      }
    }

    static std::vector<compute_t> calculate_sigma(const std::vector<compute_t>& d, const compute_t sigma_max) {
      auto sigma_bc(d);
      for (auto& x: sigma_bc) {
        x = sigma_max * std::pow(x, static_cast<compute_t>(PMLGrade));
      }
      return sigma_bc;
    }

    static std::vector<compute_t> calculate_alpha(const std::vector<compute_t>& d) {
      auto alpha_bc(d);
      for (auto& x: alpha_bc) {
        x = static_cast<compute_t>(PMLAlphaMax) * std::pow(1.0_fp - x, 1.0_fp);
      }
      return alpha_bc;
    }

    void calculate_coeffs(const std::vector<compute_t>& sigma, const std::vector<compute_t>& alpha) {
      constexpr auto coef1 = -dt / static_cast<compute_t>(constants::eps0);

      for (auto i = 0zu; i < PMLDepth; i++) {
        constexpr auto kappa_bc = 1.0_fp;
        b[i] = std::exp(coef1 * ((sigma[i] / kappa_bc) + alpha[i]));
        c[i] = (sigma[i] * (b[i] - 1.0_fp)) / (kappa_bc * (sigma[i] + (kappa_bc * alpha[i])));
      }
    }

    Array3D<compute_t> psi;
    offset_t offsets;
    coeffs_t b{};
    coeffs_t c{};
  };

  template<EMFace F, EMSide S>
  struct PeriodicData {
    explicit PeriodicData(const Array3D<compute_t>& f)
    : numInterior(getNumInterior(f.nx(), f.ny(), f.nz())),
      hi_idx(getHiIndex(f.nx(), f.ny(), f.nz())),
      offsets(getOffsets<F, S, nHalo>(f))
    {}

    static std::size_t getNumInterior(const std::size_t nx, const std::size_t ny, const std::size_t nz) {
      if constexpr (F == EMFace::X) {
        return nx - (2zu * nHalo);
      } else if constexpr (F == EMFace::Y) {
        return ny - (2zu * nHalo);
      } else {
        return nz - (2zu * nHalo);
      }
    }

    static std::size_t getHiIndex(const std::size_t nx, const std::size_t ny, const std::size_t nz) {
      if constexpr (F == EMFace::X) {
        return nx - 1 - nHalo;
      } else if constexpr (F == EMFace::Y) {
        return ny - 1 - nHalo;
      } else {
        return nz - 1 - nHalo;
      }
    }

    std::size_t numInterior;
    std::size_t hi_idx;
    std::array<std::size_t, 6> offsets;
  };
  
  template<EMFace F, EMSide S>
  struct PMLFaceBC {
    explicit PMLFaceBC() = delete;

    explicit PMLFaceBC(const auto& emdata)
    : Ex{emdata.Ex, getPMLOffsets<F, S, true>(emdata.Ex), true},
      Ey{emdata.Ey, getPMLOffsets<F, S, true>(emdata.Ey), true},
      Ez{emdata.Ez, getPMLOffsets<F, S, true>(emdata.Ez), true},
      Hx{emdata.Hx, getPMLOffsets<F, S, false>(emdata.Hx), false},
      Hy{emdata.Hy, getPMLOffsets<F, S, false>(emdata.Hy), false},
      Hz{emdata.Hz, getPMLOffsets<F, S, false>(emdata.Hz), false}
    {}

    PMLData<F, S> Ex;
    PMLData<F, S> Ey;
    PMLData<F, S> Ez;
    PMLData<F, S> Hx;
    PMLData<F, S> Hy;
    PMLData<F, S> Hz;
  };

  template<EMFace F, EMSide S>
  struct PeriodicFaceBC {
    explicit PeriodicFaceBC() = delete;

    explicit PeriodicFaceBC(const auto& emdata)
    : Ex{emdata.Ex},
      Ey{emdata.Ey},
      Ez{emdata.Ez},
      Hx{emdata.Hx},
      Hy{emdata.Hy},
      Hz{emdata.Hz}
    {}

    PeriodicData<F, S> Ex;
    PeriodicData<F, S> Ey;
    PeriodicData<F, S> Ez;
    PeriodicData<F, S> Hx;
    PeriodicData<F, S> Hy;
    PeriodicData<F, S> Hz;
  };

  template<BCType BC>
  struct BCData;
  
  template<>
  struct BCData<BCType::PML> {
    explicit BCData(const auto& emdata)
    : x0(emdata),
      y0(emdata),
      z0(emdata),
      x1(emdata),
      y1(emdata),
      z1(emdata)
    {}

    PMLFaceBC<EMFace::X, EMSide::Lo> x0;
    PMLFaceBC<EMFace::Y, EMSide::Lo> y0;
    PMLFaceBC<EMFace::Z, EMSide::Lo> z0;
    PMLFaceBC<EMFace::X, EMSide::Hi> x1;
    PMLFaceBC<EMFace::Y, EMSide::Hi> y1;
    PMLFaceBC<EMFace::Z, EMSide::Hi> z1;
  };

  template<>
  struct BCData<BCType::Periodic> {
    explicit BCData(const auto& emdata)
    : x0(emdata),
      y0(emdata),
      z0(emdata),
      x1(emdata),
      y1(emdata),
      z1(emdata)
    {}

    PeriodicFaceBC<EMFace::X, EMSide::Lo> x0;
    PeriodicFaceBC<EMFace::Y, EMSide::Lo> y0;
    PeriodicFaceBC<EMFace::Z, EMSide::Lo> z0;
    PeriodicFaceBC<EMFace::X, EMSide::Hi> x1;
    PeriodicFaceBC<EMFace::Y, EMSide::Hi> y1;
    PeriodicFaceBC<EMFace::Z, EMSide::Hi> z1;
  };
}

#endif //EM_BOUNDARIES_HPP
