#ifndef BC_FUNCTORS_HPP
#define BC_FUNCTORS_HPP

#include "diff_operators.hpp"
#include "program_params.hpp"
#include "em_params.hpp"

namespace tf::electromagnetics
{

template<typename CurlFunc, bool isLo, bool Negate>
struct PMLFunctor {
   using Curl = CurlFunc;

   #pragma omp declare simd notinbranch
   static void apply(auto& f1, const auto& f2, const auto& c1, auto& bc,
                     const auto i, const auto j, const auto k, const auto x0)
      requires (Curl::type == Derivative::DX)
   {
      auto ipml;
      if constexpr (isLo) { ipml = i; }
      else { ipml = i - x0 + Curl::Forward; }

      bc.psi(ipml, j, k) = bc.b[ipml] * bc.psi(ipml, j, k) + bc.c[ipml] * Curl::apply(f2, i, j, k);
      if constexpr (Negate) {
         f1(i, j, k) -= c1(i, j, k) * bc.psi(ipml, j, k);
      }
      else {
         f1(i, j, k) += c1(i, j, k) * bc.psi(ipml, j, k);
      }
   } // end apply

   #pragma omp declare simd notinbranch
   static void apply(auto& f1, const auto& f2, const auto& c1, auto& bc,
                     const auto i, const auto j, const auto k, const auto y0)
      requires (Curl::type == Derivative::DY)
   {
      auto jpml;
      if constexpr (isLo) { jpml = j; }
      else { jpml = j - y0 + Curl::Forward; }

      bc.psi(i, jpml, k) = bc.b[jpml] * bc.psi(i, jpml, k) + bc.c[jpml] * Curl::apply(f2, i, j, k);
      if constexpr (Negate) {
         f1(i, j, k) -= c1(i, j, k) * bc.psi(i, jpml, k);
      }
      else {
         f1(i, j, k) += c1(i, j, k) * bc.psi(i, jpml, k);
      }
   } // end apply

   #pragma omp declare simd notinbranch
   static void apply(auto& f1, const auto& f2, const auto& c1, auto& bc,
                     const auto i, const auto j, const auto k, const auto z0)
      requires (Curl::type == Derivative::DZ)
   {
      auto kpml;
      if constexpr (isLo) { kpml = k; }
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


template<Derivative D, bool Add>
struct PeriodicFunctor {
   #pragma omp declare simd notinbranch
   static void apply(auto& f, const auto& bc, const auto i, const auto j, const auto k)
   {
      auto idx1, idx2, idx3, idx4;
      if constexpr (D == Derivative::DX) {
         const auto pm = i % bc.numInterior;
         idx1 = f.get_scid(nHalo - 1 - i, j, k);
         idx2 = f.get_scid(bc.hiIndex - pm, j, k);
         idx3 = f.get_scid(bc.hiIndex + 1 + i, j, k);
         idx4 = f.get_scid(nHalo + pm, j, k);
      } else if constexpr (D == Derivative::DY) {
         const auto pm = j % bc.numInterior;
         idx1 = f.get_scid(i, nHalo - 1 - j, k);
         idx2 = f.get_scid(i, bc.hiIndex - pm, k);
         idx3 = f.get_scid(i, bc.hiIndex + 1 + j, k);
         idx4 = f.get_scid(i, nHalo + pm, k);
      } else {
         const auto pm = k % bc.numInterior;
         idx1 = f.get_scid(i, j, nHalo - 1 - k);
         idx2 = f.get_scid(i, j, bc.hiIndex - pm);
         idx3 = f.get_scid(i, j, bc.hiIndex + 1 + k);
         idx4 = f.get_scid(i, j, nHalo + pm);
      }

      // This allows for cumulative boundaries for current density
      if constexpr (Add) {
         f[idx1] += f[idx2];
         f[idx3] += f[idx4];
      } else {
         f[idx1] = f[idx2];
         f[idx3] = f[idx4];
      }
   } // end apply()
}; // end struct PeriodicFunctor


template<bool Negate>
struct BCIntegrator {
   static constexpr void apply(const auto&, const auto&, const auto&, const auto&) {}

   static void apply(auto& psi, const auto& b, const auto& c, auto& f, const auto& d, const auto& cf)
   // requires std::same_as<U, PMLData>
   {
      // #pragma omp parallel for simd collapse(3) num_threads(nThreads)
      for (auto i = 0zu; i < psi.extent(0); ++i) {
         for (auto j = 0zu; j < psi.extent(1); ++j) {
            for (auto k = 0zu; k < psi.extent(2); ++k) {
               psi[i, j, k] = b[k] * psi[i, j, k] + c[k] * (d[i, j, k + 1] - d[i, j, k]);
               f[i, j, k] += cf * psi[i, j, k]; // cf comes with -1 baked in for += or -=
            } // end for k
         } // end for j
      } // end for i
   } // end operator()

   // static void apply(auto& f1, const auto&, const auto&, const auto& bc)
   // requires std::same_as<U, PeriodicData>
   // {
   //    const auto& [x0, x1, y0, y1, z0, z1] = bc.offsets;
   //    #pragma omp parallel for simd collapse(3) num_threads(nThreads)
   //    for (auto i = x0; i < x1; ++i) {
   //       for (auto j = y0; j < y1; ++j) {
   //          for (auto k = z0; k < z1; ++k) {
   //             // todo: Need to make this work with J periodic, aka add = true
   //             PeriodicFunctor<Curl::type, false>::apply(f1, bc, i, j, k);
   //          }
   //       }
   //    }
   // }
}; // end struct BCIntegrator

} // end namespace tf::electromagnetics


#endif //BC_FUNCTORS_HPP
