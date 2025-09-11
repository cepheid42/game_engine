#ifndef BC_FUNCTORS_HPP
#define BC_FUNCTORS_HPP

#include "program_params.hpp"
#include "em_params.hpp"
#include "diff_operators.hpp"
#include "traits.hpp"

#include <concepts>

namespace tf::electromagnetics {


template<typename T> concept is_periodic = std::derived_from<T, periodic_t>;
template<typename T> concept is_pml = std::derived_from<T, pml_t>;

template<typename CurlFunc, bool isLo, bool Negate>
struct PMLFunctor : pml_t {
   using Curl                   = CurlFunc;
   // static constexpr auto ishi   = isHi;
   // static constexpr auto negate = Negate;

   #pragma omp declare simd notinbranch
   static void apply(auto& f1, const auto& f2, const auto& c1, auto& bc,
                     const std::size_t i, const std::size_t j, const std::size_t k, const std::size_t x0)
      requires (Curl::type == Derivative::DX)
   {
      std::size_t ipml;
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
                     const std::size_t i, const std::size_t j, const std::size_t k, const std::size_t y0)
      requires (Curl::type == Derivative::DY)
   {
      std::size_t jpml;
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
                     const std::size_t i, const std::size_t j, const std::size_t k, const std::size_t z0)
      requires (Curl::type == Derivative::DZ)
   {
      std::size_t kpml;
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


template<EMFace F, bool Add>
struct PeriodicFunctor : periodic_t {
   // static constexpr EMFace face = F;
   // static constexpr bool add = Add;

   static void apply(auto& f, const auto& bc, const std::size_t i, const std::size_t j, const std::size_t k)
   {
      std::size_t idx1, idx2, idx3, idx4;
      if constexpr (F == EMFace::X) {
         const auto pm = i % bc.numInterior;
         idx1 = f.get_scid(NHalo - 1 - i, j, k);
         idx2 = f.get_scid(bc.hiIndex - pm, j, k);
         idx3 = f.get_scid(bc.hiIndex + 1 + i, j, k);
         idx4 = f.get_scid(NHalo + pm, j, k);
      } else if constexpr (F == EMFace::Y) {
         const auto pm = j % bc.numInterior;
         idx1 = f.get_scid(i, NHalo - 1 - j, k);
         idx2 = f.get_scid(i, bc.hiIndex - pm, k);
         idx3 = f.get_scid(i, bc.hiIndex + 1 + j, k);
         idx4 = f.get_scid(i, NHalo + pm, k);
      } else {
         const auto pm = k % bc.numInterior;
         idx1 = f.get_scid(i, j, NHalo - 1 - k);
         idx2 = f.get_scid(i, j, bc.hiIndex - pm);
         idx3 = f.get_scid(i, j, bc.hiIndex + 1 + k);
         idx4 = f.get_scid(i, j, NHalo + pm);
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


template<typename U>
struct BCIntegrator {
   static constexpr void apply(const auto&, const auto&, const auto&, const auto&) {}

   static void apply(auto& f1, const auto& f2, const auto& c1, auto& bc)
   requires is_pml<U>
   {
      static constexpr Derivative CurlType = U::Curl::type;
      const auto& [x0, x1, y0, y1, z0, z1] = bc.offsets;
   
      std::size_t pml_offset;
      if constexpr      (CurlType == Derivative::DX) { pml_offset = x0; }
      else if constexpr (CurlType == Derivative::DY) { pml_offset = y0; }
      else                                           { pml_offset = z0; }
   
      #pragma omp parallel for simd collapse(3) num_threads(nThreads)
      for (std::size_t i = x0; i < x1; ++i) {
         for (std::size_t j = y0; j < y1; ++j) {
            for (std::size_t k = z0; k < z1; ++k) {
               U::apply(f1, f2, c1, bc, i, j, k, pml_offset);
            } // end for k
         } // end for j
      } // end for i
   } // end operator()
   
   static void apply(auto& f1, const auto&, const auto&, const auto& bc)
   requires is_periodic<U>
   {
      const auto& [x0, x1, y0, y1, z0, z1] = bc.offsets;
      #pragma omp parallel for simd collapse(3) num_threads(nThreads)
      for (std::size_t i = x0; i < x1; ++i) {
         for (std::size_t j = y0; j < y1; ++j) {
            for (std::size_t k = z0; k < z1; ++k) {
               U::apply(f1, bc, i, j, k);
            }
         }
      }
   }
}; // end struct BCIntegrator

template<EMFace F, EMSide S>
struct PeriodicImpl;

template<EMFace F>
struct PeriodicImpl<F, EMSide::Hi> {
   using Ex = BCIntegrator<null_t>;
   using Ey = BCIntegrator<null_t>;
   using Ez = BCIntegrator<null_t>;
   using Hx = BCIntegrator<null_t>;
   using Hy = BCIntegrator<null_t>;
   using Hz = BCIntegrator<null_t>;
   using Jx = BCIntegrator<null_t>;
   using Jy = BCIntegrator<null_t>;
   using Jz = BCIntegrator<null_t>;
};

template<>
struct PeriodicImpl<EMFace::X, EMSide::Lo> {
   using Ex = BCIntegrator<null_t>;
   using Ey = BCIntegrator<PeriodicFunctor<EMFace::X, false>>;
   using Ez = BCIntegrator<PeriodicFunctor<EMFace::X, false>>;
   using Hx = BCIntegrator<null_t>;
   using Hy = BCIntegrator<null_t>;
   using Hz = BCIntegrator<null_t>;
   using Jx = BCIntegrator<null_t>;
   using Jy = BCIntegrator<null_t>;
   using Jz = BCIntegrator<null_t>;
};

template<>
struct PeriodicImpl<EMFace::Y, EMSide::Lo> {
   using Ex = BCIntegrator<PeriodicFunctor<EMFace::Y, false>>;
   using Ey = BCIntegrator<null_t>;
   using Ez = BCIntegrator<PeriodicFunctor<EMFace::Y, false>>;
   using Hx = BCIntegrator<null_t>;
   using Hy = BCIntegrator<null_t>;
   using Hz = BCIntegrator<null_t>;
   using Jx = BCIntegrator<null_t>;
   using Jy = BCIntegrator<null_t>;
   using Jz = BCIntegrator<null_t>;
};

template<>
struct PeriodicImpl<EMFace::Z, EMSide::Lo> {
   using Ex = BCIntegrator<PeriodicFunctor<EMFace::Z, false>>;
   using Ey = BCIntegrator<PeriodicFunctor<EMFace::Z, false>>;
   using Ez = BCIntegrator<null_t>;
   using Hx = BCIntegrator<null_t>;
   using Hy = BCIntegrator<null_t>;
   using Hz = BCIntegrator<null_t>;
   using Jx = BCIntegrator<null_t>;
   using Jy = BCIntegrator<null_t>;
   using Jz = BCIntegrator<null_t>;
};

template<EMFace F, EMSide S>
struct PMLImpl;

template<EMSide S>
struct PMLImpl<EMFace::X, S> {
   static constexpr bool isLo = S == EMSide::Lo;
   using Ex = BCIntegrator<null_t>;
   using Ey = BCIntegrator<PMLFunctor<backward_dx, isLo, true>>;
   using Ez = BCIntegrator<PMLFunctor<backward_dx, isLo, false>>;
   using Hx = BCIntegrator<null_t>;
   using Hy = BCIntegrator<PMLFunctor<forward_dx, isLo, false>>;
   using Hz = BCIntegrator<PMLFunctor<forward_dx, isLo, true>>;
   using Jx = BCIntegrator<null_t>;
   using Jy = BCIntegrator<null_t>;
   using Jz = BCIntegrator<null_t>;
};

template<EMSide S>
struct PMLImpl<EMFace::Y, S> {
   static constexpr bool isLo = S == EMSide::Lo;
   using Ex = BCIntegrator<PMLFunctor<backward_dy, isLo, false>>;
   using Ey = BCIntegrator<null_t>;
   using Ez = BCIntegrator<PMLFunctor<backward_dy, isLo, true>>;
   using Hx = BCIntegrator<PMLFunctor<forward_dy, isLo, true>>;
   using Hy = BCIntegrator<null_t>;
   using Hz = BCIntegrator<PMLFunctor<forward_dy, isLo, false>>;
   using Jx = BCIntegrator<null_t>;
   using Jy = BCIntegrator<null_t>;
   using Jz = BCIntegrator<null_t>;
};

template<EMSide S>
struct PMLImpl<EMFace::Z, S> {
   static constexpr bool isLo = S == EMSide::Lo;
   using Ex = BCIntegrator<PMLFunctor<backward_dz, isLo, true>>;
   using Ey = BCIntegrator<PMLFunctor<backward_dz, isLo, false>>;
   using Ez = BCIntegrator<null_t>;
   using Hx = BCIntegrator<PMLFunctor<forward_dz, isLo, false>>;
   using Hy = BCIntegrator<PMLFunctor<forward_dz, isLo, true>>;
   using Hz = BCIntegrator<null_t>;
   using Jx = BCIntegrator<null_t>;
   using Jy = BCIntegrator<null_t>;
   using Jz = BCIntegrator<null_t>;
};

struct ReflectingBoundary {
    using Ex = BCIntegrator<null_t>;
    using Ey = BCIntegrator<null_t>;
    using Ez = BCIntegrator<null_t>;
    using Hx = BCIntegrator<null_t>;
    using Hy = BCIntegrator<null_t>;
    using Hz = BCIntegrator<null_t>;
    using Jx = BCIntegrator<null_t>;
    using Jy = BCIntegrator<null_t>;
    using Jz = BCIntegrator<null_t>;
};

template<EMFace F, EMSide S>
using PeriodicBoundary = PeriodicImpl<F, S>;

template<EMFace F, EMSide S>
using PMLBoundary = PMLImpl<F, S>;

} // end namespace tf::electromagnetics


#endif //BC_FUNCTORS_HPP
