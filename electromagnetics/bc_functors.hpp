#ifndef BC_FUNCTORS_HPP
#define BC_FUNCTORS_HPP

#include "diff_operators.hpp"

#include <concepts>

namespace tf::electromagnetics {
template<typename CurlFunc, bool Hi, bool Negate>
struct PMLFunctor {
   using Curl                   = CurlFunc;
   static constexpr auto hi     = Hi;
   static constexpr auto negate = Negate;

   static void apply(auto& f1, const auto& f2, const auto& c1, auto& bc, const std::size_t i, const std::size_t j, const std::size_t k, const std::size_t x0)
      requires (Curl::type == Derivative::DX)
   {
      std::size_t ipml;
      if constexpr (!Hi) { ipml = i; }
      else { ipml = i - x0 + Curl::Forward; }

      bc.psi(ipml, j, k) = bc.b[ipml] * bc.psi(ipml, j, k) + bc.c[ipml] * Curl::apply(f2, i, j, k);
      if constexpr (Negate) {
         f1(i, j, k) -= c1(i, j, k) * bc.psi(ipml, j, k);
      }
      else {
         f1(i, j, k) += c1(i, j, k) * bc.psi(ipml, j, k);
      }
   } // end apply

   static void apply(auto& f1, const auto& f2, const auto& c1, auto& bc, const std::size_t i, const std::size_t j, const std::size_t k, const std::size_t y0)
      requires (Curl::type == Derivative::DY)
   {
      std::size_t jpml;
      if constexpr (!Hi) { jpml = j; }
      else { jpml = j - y0 + Curl::Forward; }

      bc.psi(i, jpml, k) = bc.b[jpml] * bc.psi(i, jpml, k) + bc.c[jpml] * Curl::apply(f2, i, j, k);
      if constexpr (Negate) {
         f1(i, j, k) -= c1(i, j, k) * bc.psi(i, jpml, k);
      }
      else {
         f1(i, j, k) += c1(i, j, k) * bc.psi(i, jpml, k);
      }
   } // end apply

   static void apply(auto& f1, const auto& f2, const auto& c1, auto& bc, const std::size_t i, const std::size_t j, const std::size_t k, const std::size_t z0)
      requires (Curl::type == Derivative::DZ)
   {
      std::size_t kpml;
      if constexpr (!Hi) { kpml = k; }
      else { kpml = k - z0 + Curl::Forward; }

      bc.psi(i, j, kpml) = bc.b[kpml] * bc.psi(i, j, kpml) + bc.c[kpml] * Curl::apply(f2, i, j, k);
      if constexpr (Negate) {
         f1(i, j, k) -= c1(i, j, k) * bc.psi(i, j, kpml);
      }
      else {
         f1(i, j, k) += c1(i, j, k) * bc.psi(i, j, kpml);
      }
   } // end apply
};   // end struct Pml3DFunctor


template<typename UpdateFunc>
struct BCIntegrator {
   static void operator()(auto& f1, const auto& f2, const auto& c1, auto& bc)
      requires std::same_as<UpdateFunc, PMLFunctor<typename UpdateFunc::Curl, UpdateFunc::hi, UpdateFunc::negate>>
   {
      const auto& [x0, x1, y0, y1, z0, z1] = bc.offsets;
      std::size_t pml_offset;
      if constexpr (UpdateFunc::Curl::type == Derivative::DX) { pml_offset = x0; }
      else if constexpr (UpdateFunc::Curl::type == Derivative::DY) { pml_offset = y0; }
      else { pml_offset = z0; }

      #pragma omp parallel for simd collapse(3) num_threads(nThreads)
      for (std::size_t i = x0; i < x1; ++i) {
         for (std::size_t j = y0; j < y1; ++j) {
            for (std::size_t k = z0; k < z1; ++k) {
               UpdateFunc::apply(f1, f2, c1, bc, i, j, k, pml_offset);
            } // end for k
         } // end for j
      } // end for i
   } // end operator()
}; // end struct BCIntegrator

template<>
struct BCIntegrator<void> {
   static constexpr void operator()() {}
   static constexpr void operator()(const auto&, const auto&, const auto&, const auto&) {}
};
} // end namespace tf::electromagnetics


#endif //BC_FUNCTORS_HPP
