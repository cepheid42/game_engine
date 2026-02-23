#ifndef EM_BOUNDARIES_HPP
#define EM_BOUNDARIES_HPP

#include "program_params.hpp"
#include "array.hpp"
#include "constants.hpp"
#include "math_utils.hpp"

#include <vector>
#include <algorithm>
#include <array>
#include <cassert>

namespace tf::electromagnetics {
template<EMSide S, std::size_t Depth>
constexpr std::array<std::size_t, 6> get_x_offsets(const std::size_t nx, const std::size_t ny, const std::size_t nz) {
   if constexpr (S == EMSide::Lo) { return {0, Depth, 0, ny, 0, nz}; }
   else                           { return {nx - Depth, nx, 0, ny, 0, nz}; }
}

template<EMSide S, std::size_t Depth>
constexpr std::array<std::size_t, 6> get_y_offsets(const std::size_t nx, const std::size_t ny, const std::size_t nz) {
   if constexpr (S == EMSide::Lo) { return {0, nx, 0, Depth, 0, nz}; }
   else                           { return {0, nx, ny - Depth, ny, 0, nz}; }
}

template<EMSide S, std::size_t Depth>
constexpr std::array<std::size_t, 6> get_z_offsets(const std::size_t nx, const std::size_t ny, const std::size_t nz) {
   if constexpr (S == EMSide::Lo) { return {0, nx, 0, ny, 0, Depth}; }
   else                           { return {0, nx, 0, ny, nz - Depth, nz};  }
}

template<EMFace F, EMSide S, std::size_t Depth>
constexpr auto getOffsets(const auto& f) {
   if constexpr (F == EMFace::X) {
      return get_x_offsets<S, Depth>(f.nx(), f.ny(), f.nz());
   }
   else if constexpr (F == EMFace::Y) {
      return get_y_offsets<S, Depth>(f.nx(), f.ny(), f.nz());
   }
   else {
      return get_z_offsets<S, Depth>(f.nx(), f.ny(), f.nz());
   }
} // end getOffsets()

template<EMFace F, EMSide S, bool isE>
auto getPMLOffsets(const auto& f) {
   auto           result = getOffsets<F, S, PMLDepth>(f);
   constexpr auto lo     = static_cast<std::size_t>(isE == (S == EMSide::Lo)); // xnor
   constexpr auto hi     = static_cast<std::size_t>(isE != (S == EMSide::Lo)); // xor

   std::size_t lo_idx, hi_idx;
   if constexpr (F == EMFace::X) {
      lo_idx = 0; // x0
      hi_idx = 1; // x1
   }
   else if constexpr (F == EMFace::Y) {
      lo_idx = 2; // y0
      hi_idx = 3; // y1
   }
   else {
      lo_idx = 4; // z0
      hi_idx = 5; // z1
   }

   result[lo_idx] += lo;
   result[hi_idx] -= hi;

   return result;
} // end getPMLOffsets()


template<EMFace F, EMSide S, bool isE>
struct PMLData {
   using offset_t = std::array<std::size_t, 6>;
   using coeffs_t = std::array<double, PMLDepth>;

   explicit PMLData() = delete;

   PMLData(const std::size_t nx_, const std::size_t ny_, const std::size_t nz_, const auto delta, const offset_t& offsets_)
   : psi(nx_, ny_, nz_),
     offsets(offsets_)
   { init_coefficients(delta); }

   PMLData(const Array3D<double>& f, const offset_t& offset) requires (F == EMFace::X)
   : PMLData(PMLDepth, f.ny(), f.nz(), dx, offset) {}

   PMLData(const Array3D<double>& f, const offset_t& offset) requires (F == EMFace::Y)
   : PMLData(f.nx(), PMLDepth, f.nz(), dy, offset) {}

   PMLData(const Array3D<double>& f, const offset_t& offset) requires (F == EMFace::Z)
   : PMLData(f.nx(), f.ny(), PMLDepth, dz, offset) {}

   void init_coefficients(const auto delta) {
      std::vector<double> d = math::linspace(1.0, 0.0, PMLDepth, false);

      if (!isE) {
         constexpr auto hstep = 1.0 / (2.0 * static_cast<double>(PMLDepth));
         for (auto& x: d) { x -= hstep; }
      }

      constexpr auto eta0      = constants::eta0;
      constexpr auto pml_grade = static_cast<double>(PMLGrade);
      const auto sigma_max = 0.8 * (pml_grade + 1.0) / (delta * eta0);

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
      for (auto& x: sigma_bc) { x = sigma_max * std::pow(x, PMLGrade); }
      return sigma_bc;
   }

   static std::vector<double> calculate_alpha(const std::vector<double>& d) {
      static constexpr auto m = 1.0; // alpha exponent
      auto alpha_bc(d);
      for (auto& x: alpha_bc) { x = PMLAlphaMax * std::pow(1.0 - x, m); }
      return alpha_bc;
   }

    void calculate_coeffs(const std::vector<double>& sigma, const std::vector<double>& alpha) {
       static constexpr auto kappa_bc = 1.0;
       constexpr auto coef1 = -dt_em / constants::eps0;
       for (auto i = 0lu; i < PMLDepth; i++) {
        b[i] = std::exp(coef1 * ((sigma[i] / kappa_bc) + alpha[i]));
        c[i] = (sigma[i] * (b[i] - 1.0)) / (kappa_bc * (sigma[i] + (kappa_bc * alpha[i])));
      }
   }

   Array3D<double> psi;
   offset_t        offsets;
   coeffs_t        b{};
   coeffs_t        c{};
}; // end struct PMLData


template<EMFace F, EMSide S>
struct PeriodicData {
   using offset_t = std::array<std::size_t, 6>;

   explicit PeriodicData(const auto& f)
   : numInterior{get_num_interior(f)},
     hiIndex{get_hi_index(f)},
     offsets{getOffsets<F, S, nHalo>(f)}
   {}

   static std::size_t get_num_interior(const auto& f) {
      if      constexpr (F == EMFace::X) { return f.nx() - (2 * nHalo); }
      else if constexpr (F == EMFace::Y) { return f.ny() - (2 * nHalo); }
      else                               { return f.nz() - (2 * nHalo); }
   }

   static std::size_t get_hi_index(const auto& f) {
      if      constexpr (F == EMFace::X) { return f.nx() - 1 - nHalo; }
      else if constexpr (F == EMFace::Y) { return f.ny() - 1 - nHalo; }
      else                               { return f.nz() - 1 - nHalo; }
   }
   
   std::size_t numInterior;
   std::size_t hiIndex;
   offset_t offsets;
};

struct EmptyData {};

template<EMFace F, EMSide S>
struct PMLFaceBC {
   explicit PMLFaceBC(const auto& emdata)
   : Ex{emdata.Ex, getPMLOffsets<F, S, true>(emdata.Ex)},
     Ey{emdata.Ey, getPMLOffsets<F, S, true>(emdata.Ey)},
     Ez{emdata.Ez, getPMLOffsets<F, S, true>(emdata.Ez)},
     Hx{emdata.Hx, getPMLOffsets<F, S, false>(emdata.Hx)},
     Hy{emdata.Hy, getPMLOffsets<F, S, false>(emdata.Hy)},
     Hz{emdata.Hz, getPMLOffsets<F, S, false>(emdata.Hz)}
   {}

   PMLData<F, S, true> Ex;
   PMLData<F, S, true> Ey;
   PMLData<F, S, true> Ez;
   PMLData<F, S, false> Hx;
   PMLData<F, S, false> Hy;
   PMLData<F, S, false> Hz;
   EmptyData Jx{};
   EmptyData Jy{};
   EmptyData Jz{};
}; // end struct PMLFaceBC


template<EMFace F, EMSide S>
struct PeriodicFaceBC {
   explicit PeriodicFaceBC(const auto& emdata)
   : Ex{emdata.Ex}, Ey{emdata.Ey}, Ez{emdata.Ez},
     Hx{emdata.Hx}, Hy{emdata.Hy}, Hz{emdata.Hz},
     Jx{emdata.Jx}, Jy{emdata.Jy}, Jz{emdata.Jz}
   {}

   PeriodicData<F, S> Ex;
   PeriodicData<F, S> Ey;
   PeriodicData<F, S> Ez;
   PeriodicData<F, S> Hx;
   PeriodicData<F, S> Hy;
   PeriodicData<F, S> Hz;
   PeriodicData<F, S> Jx;
   PeriodicData<F, S> Jy;
   PeriodicData<F, S> Jz;
}; // end struct PeriodicFaceBC

template<EMFace, EMSide>
struct ReflectingFaceBC {
   explicit ReflectingFaceBC(const auto&) {}
   EmptyData Ex{};
   EmptyData Ey{};
   EmptyData Ez{};
   EmptyData Hx{};
   EmptyData Hy{};
   EmptyData Hz{};
   EmptyData Jx{};
   EmptyData Jy{};
   EmptyData Jz{};
};

template<
   typename X0BC, typename X1BC,
   typename Y0BC, typename Y1BC,
   typename Z0BC, typename Z1BC
>
struct BCData {
   explicit BCData(const auto& emdata)
   : x0(emdata), y0(emdata), z0(emdata),
     x1(emdata), y1(emdata), z1(emdata)
   {}

   X0BC x0;
   Y0BC y0;
   Z0BC z0;
   X1BC x1;
   Y1BC y1;
   Z1BC z1;
}; // end struct BCData
} // end namespace tf::electromagnetics

#endif //EM_BOUNDARIES_HPP
