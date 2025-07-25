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

   #pragma omp declare simd notinbranch
   static void apply(auto& f1, const auto& f2, const auto& c1, auto& bc,
                     const std::size_t i, const std::size_t j, const std::size_t k, const std::size_t x0)
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

   #pragma omp declare simd notinbranch
   static void apply(auto& f1, const auto& f2, const auto& c1, auto& bc,
                     const std::size_t i, const std::size_t j, const std::size_t k, const std::size_t y0)
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

   #pragma omp declare simd notinbranch
   static void apply(auto& f1, const auto& f2, const auto& c1, auto& bc,
                     const std::size_t i, const std::size_t j, const std::size_t k, const std::size_t z0)
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

// template<EMFace F, EMSide S>
// struct PeriodicFunctor {
//    static constexpr EMFace face = F;
//    static constexpr EMSide side = S;
//
//    // static void updateH() { DBG("Periodic3D::updateE()"); }
//    static void updateE(auto&, auto&) {} // Only lo-side periodic is used
//
//    static void updateE(auto& f, const auto& bc, const std::size_t i, const std::size_t j, const std::size_t k)
//       requires (S == EMSide::Lo)
//    {
//       const auto& numInterior = bc.numInterior;
//       const auto& hiIndex = bc.hiIndex;
//       if constexpr (F == EMFace::X) {
//          const auto pm = i % numInterior;
//          f(bc.depth - 1 - i, j, k) = f(hiIndex - pm, j, k);
//          f(hiIndex + 1 + i, j, k) = f(bc.depth + pm, j, k);
//       } else if constexpr (F == EMFace::Y) {
//          const auto pm = j % numInterior;
//          f(i, bc.depth - 1 - j, k) = f(i, hiIndex - pm, k);
//          f(i, hiIndex + 1 + j, k) = f(i, bc.depth + pm, k);
//       } else {
//          const auto pm = k % numInterior;
//          f(i, j, bc.depth - 1 - k) = f(i, j, hiIndex - pm);
//          f(i, j, hiIndex + 1 + k) = f(i, j, bc.depth + pm);
//       }
//    } // end updateH
// };

template<EMFace F, bool Add=false>
struct PeriodicFunctor {
   static constexpr auto nHalo = 2zu;

   static void apply(auto& f, const auto& bc, const std::size_t i, const std::size_t j, const std::size_t k)
   {
      std::size_t idx1, idx2, idx3, idx4;
      if constexpr (F == EMFace::X) {
         const auto pm = i % bc.numInterior;
         idx1 = f.get_scid(nHalo - 1 - i, j, k);
         idx2 = f.get_scid(bc.hiIndex - pm, j, k);
         idx3 = f.get_scid(bc.hiIndex + 1 + i, j, k);
         idx4 = f.get_scid(nHalo + pm, j, k);
      } else if constexpr (F == EMFace::Y) {
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

      if constexpr (Add) {
         f[idx1] += f[idx2];
         f[idx3] += f[idx4];
      } else {
         f[idx1] = f[idx2];
         f[idx3] = f[idx4];
      }
   } // end apply()
}; // end struct PeriodicFunctor

template<typename UpdateFunc>
struct BCIntegrator {
   static void operator()(auto& f1, const auto& bc) {
      const auto& [x0, x1, y0, y1, z0, z1] = bc.offsets;
      // #pragma omp parallel for simd collapse(3) num_threads(nThreads)
      for (std::size_t i = x0; i < x1; ++i) {
         for (std::size_t j = y0; j < y1; ++j) {
            for (std::size_t k = z0; k < z1; ++k) {
               UpdateFunc::apply(f1, bc, i, j, k);
            }
         }
      }
   }

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
   static constexpr void operator()(const auto&, const auto&) {}
   static constexpr void operator()(const auto&, const auto&, const auto&, const auto&) {}
};
} // end namespace tf::electromagnetics


#endif //BC_FUNCTORS_HPP
