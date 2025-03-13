#ifndef BC_FUNCTORS_HPP
#define BC_FUNCTORS_HPP

#include "diff_operators.hpp"

#include <array>

namespace tf::electromagnetics {
  template<typename Curl, bool Hi, bool Negate, bool Forward>
  struct PMLFunctor {
    static void apply(auto& f1, const auto& f2, const auto& c1, auto& psi, const auto& b, const auto& c, const std::size_t i, const std::size_t j, const std::size_t k, const std::size_t x0)
    requires (Curl::type == Derivative::DX)
    {
      std::size_t ipml;
      if constexpr (!Hi) { ipml = i; }
      else { ipml = i - x0 + Forward; }

      psi(ipml, j, k) = b[ipml] * psi(ipml, j, k) + c[ipml] * Curl::apply(f2, i, j, k);
      if constexpr (Negate) {
        f1(i, j, k) -= c1(i, j, k) * psi(ipml, j, k);
      } else {
        f1(i, j, k) += c1(i, j, k) * psi(ipml, j, k);
      }
    } // end apply

    static void apply(auto& f1, const auto& f2, const auto& c1, auto& psi, const auto& b, const auto& c, const std::size_t i, const std::size_t j, const std::size_t k, const std::size_t y0)
    requires (Curl::type == Derivative::DY)
    {
      std::size_t jpml;
      if constexpr (!Hi) { jpml = j; }
      else { jpml = j - y0 + Forward; }

      psi(i, jpml, k) = b[jpml] * psi(i, jpml, k) + c[jpml] * Curl::apply(f2, i, j, k);
      if constexpr (Negate) {
        f1(i, j, k) -= c1(i, j, k) * psi(i, jpml, k);
      } else {
        f1(i, j, k) += c1(i, j, k) * psi(i, jpml, k);
      }
    } // end apply

    static void apply(auto& f1, const auto& f2, const auto& c1, auto& psi, const auto& b, const auto& c, const std::size_t i, const std::size_t j, const std::size_t k, const std::size_t z0)
    requires (Curl::type == Derivative::DZ)
    {
      std::size_t kpml;
      if constexpr (!Hi) { kpml = k; }
      else { kpml = k - z0 + Forward; }

      psi(i, j, kpml) = b[kpml] * psi(i, j, kpml) + c[kpml] * Curl::apply(f2, i, j, k);
      if constexpr (Negate) {
        f1(i, j, k) -= c1(i, j, k) * psi(i, j, kpml);
      } else {
        f1(i, j, k) += c1(i, j, k) * psi(i, j, kpml);
      }
    } // end apply
  }; // end struct Pml3DFunctor

  template<typename T, typename UpdateFunc>
  struct BCIntegrator {
    using offset_t = std::array<std::size_t, 6>;
    static constexpr auto direction = UpdateFunc::Curl::type;

    static void operator()(auto& f1, const auto& f2, const auto& c1, auto& psi, const auto& b, const auto& c, const offset_t& offsets)
    {
      const auto& [x0, x1, y0, y1, z0, z1] = offsets;
      std::size_t pml_offset;
      if constexpr      (direction == Derivative::DX) { pml_offset = x0; }
      else if constexpr (direction == Derivative::DY) { pml_offset = y0; }
      else                                            { pml_offset = z0; }

      for (std::size_t i = x0; i < x1; ++i) {
        for (std::size_t j = y0; j < y1; ++j) {
          for (std::size_t k = z0; k < z1; ++k) {
            UpdateFunc::apply(f1, f2, c1, psi, b, c, i, j, k, pml_offset);
          } // end for k
        } // end for j
      } // end for i
    } // end operator()
  }; // end struct BCIntegrator
} // end namespace tf::electromagnetics

#endif //BC_FUNCTORS_HPP
