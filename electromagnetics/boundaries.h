//
// Created by cepheid on 10/17/24.
//

#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include "electromagnetics.param"
#include "em_traits.h"

namespace tf::electromagnetics::boundaries
{
  template<typename UpdateFunc>
  struct Boundary {
    static void updateH(auto& bc, auto& f1, auto& f2, const auto& c1) {
      if constexpr (is_periodic<UpdateFunc>) {
        UpdateFunc::updateE(bc, f1);
      } else if constexpr (is_pml<UpdateFunc>) {
        UpdateFunc::updateE(bc, f1, f2, c1);
      }
      // Reflecting do not update the H-field components
    }

    static void updateE(auto& bc, auto& f1, auto& f2, const auto& c1) {
      if constexpr (is_pml<UpdateFunc>) {
        UpdateFunc::updateH(bc, f1, f2, c1);
      }
      // Reflecting/Periodic does not update E field components
    }
  };

  struct ReflectingBCUpdate {};

  template<EMFace F, EMSide S>
  struct Periodic3DUpdate : periodic_t<F, S> {
    static void updateE(auto&, auto&) {} // Only lo-side periodic is used
    static void updateE(auto& bc, auto& f) requires (S == EMSide::Lo) {
      const auto& os = bc.offsets;
      const auto& numInterior = bc.numInterior;
      const auto& hi_idx = bc.hi_idx;

      for (size_t i = os.x0; i < os.x1; ++i) {
        for (size_t j = os.y0; j < os.y1; ++j) {
          for (size_t k = os.z0; k < os.z1; ++k) {
            if constexpr (F == EMFace::X) {
              const auto pm = i % numInterior;

              f(nHalo - 1 - i, j, k) = f(hi_idx - pm, j, k);
              f(hi_idx + 1 + i, j, k) = f(nHalo + pm, j, k);
            } else if constexpr (F == EMFace::Y) {
              const auto pm = j % numInterior;

              f(i, nHalo - 1 - j, k) = f(i, hi_idx - pm, k);
              f(i, hi_idx + 1 + j, k) = f(i, nHalo + pm, k);
            } else {
              const auto pm = k % numInterior;

              f(i, j, nHalo - 1 - k) = f(i, j, hi_idx - pm);
              f(i, j, hi_idx + 1 + k) = f(i, j, nHalo + pm);
            }
          } // end k-loop
        } // end j-loop
      } // end i-loop
    } // end updateH
  };

  template<EMFace F, EMSide S, Derivative D1, bool Negate, bool Forward>
  struct Pml3DFunctor {
    using Curl = Diff<D1, Forward>;

    static void apply(auto& f1, const auto& f2, const auto& c1, auto& psi, const auto& b, const auto& c, const size_t i, const size_t j, const size_t k, const size_t x0)
    requires (F == EMFace::X)
    {
      size_t ipml;
      if constexpr (S == EMSide::Lo) { ipml = i; }
      else { ipml = i - x0 + Forward; }

      psi(ipml, j, k) = b[ipml] * psi(ipml, j, k) + c[ipml] * Curl::apply(f2, i, j, k);
      if constexpr (Negate) {
        f1(i, j, k) -= c1(i, j, k) * psi(ipml, j, k);
      } else {
        f1(i, j, k) += c1(i, j, k) * psi(ipml, j, k);
      }
    }

    static void apply(auto& f1, const auto& f2, const auto& c1, auto& psi, const auto& b, const auto& c, const size_t i, const size_t j, const size_t k, const size_t y0)
    requires (F == EMFace::Y)
    {
      size_t jpml;
      if constexpr (S == EMSide::Lo) { jpml = j; }
      else { jpml = j - y0 + Forward; }

      psi(i, jpml, k) = b[jpml] * psi(i, jpml, k) + c[jpml] * Curl::apply(f2, i, j, k);
      if constexpr (Negate) {
        f1(i, j, k) -= c1(i, j, k) * psi(i, jpml, k);
      } else {
        f1(i, j, k) += c1(i, j, k) * psi(i, jpml, k);
      }
    }

    static void apply(auto& f1, const auto& f2, const auto& c1, auto& psi, const auto& b, const auto& c, const size_t i, const size_t j, const size_t k, const size_t z0)
    requires (F == EMFace::Z)
    {
      size_t kpml;
      if constexpr (S == EMSide::Lo) { kpml = k; }
      else { kpml = k - z0 + Forward; }

      psi(i, j, kpml) = b[kpml] * psi(i, j, kpml) + c[kpml] * Curl::apply(f2, i, j, k);
      if constexpr (Negate) {
        f1(i, j, k) -= c1(i, j, k) * psi(i, j, kpml);
      } else {
        f1(i, j, k) += c1(i, j, k) * psi(i, j, kpml);
      }
    }
  };

  template<EMFace F, EMSide S, Derivative D1, bool Negate, bool isE>
  struct BCIntegrator3D {
    static auto apply(auto& f1, const auto& f2, const auto& c1, auto& bc) {
      size_t pml_offset;
      if constexpr      (F == EMFace::X) { pml_offset = bc.offsets.x0; }
      else if constexpr (F == EMFace::Y) { pml_offset = bc.offsets.y0; }
      else                               { pml_offset = bc.offsets.z0; }

#pragma omp parallel for simd collapse(3) num_threads(NTHREADS)
      for (size_t i = bc.offsets.x0; i < bc.offsets.x1; ++i) {
        for (size_t j = bc.offsets.y0; j < bc.offsets.y1; ++j) {
          for (size_t k = bc.offsets.z0; k < bc.offsets.z1; ++k) {
            Pml3DFunctor<F, S, D1, Negate, isE>::apply(f1, f2, c1, bc.psi, bc.b, bc.c, i, j, k, pml_offset);
          } // end k-loop
        } // end j-loop
      } // end i-loop
    }
  };

  template<EMFace F, EMSide S, Derivative D1, bool Negate, bool isE>
  struct Pml3DUpdate : pml_t<F, S> {
    static void updateE(auto& bc, auto& f1, const auto& f2, const auto& c1) {
      BCIntegrator3D<F, S, D1, Negate, isE>::apply(f1, f2, c1, bc);
    }

    static void updateH(auto& bc, auto& f1, const auto& f2, const auto& c1) {
      BCIntegrator3D<F, S, D1, Negate, isE>::apply(f1, f2, c1, bc);
    }
  };
} // end namespace tf::electromagnetics::boundaries

#endif //BOUNDARIES_H
