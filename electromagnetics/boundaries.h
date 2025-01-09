//
// Created by cepheid on 10/17/24.
//

#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include "em_traits.h"

namespace tf::electromagnetics::boundaries
{
  template<typename UpdateFunc>
  struct Boundary {
    static void updateE(auto& bc, auto& f1, auto& f2, const auto& c1) {
      if constexpr (is_periodic<UpdateFunc>) {
        UpdateFunc::updateE(bc, f1);
      } else if constexpr (is_pml<UpdateFunc>) {
        UpdateFunc::updateE(bc, f1, f2, c1);
      }
      // Reflecting do not update the E-field components
    }

    static void updateH(auto& bc, auto& f1, auto& f2, const auto& c1) {
      // if constexpr (is_periodic<UpdateFunc>) {
      //   UpdateFunc::updateH(bc, f1);
      // } else
      if constexpr (is_pml<UpdateFunc>) {
        UpdateFunc::updateH(bc, f1, f2, c1);
      }
      // Reflecting does not update H field components
    }
  };

  struct ReflectingBCUpdate {};

  template<EMFace F, EMSide S>
  struct Periodic3DUpdate : periodic_t<F, S> {
    // static void updateH() { DBG("Periodic3D::updateE()"); }
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

              f(bc.depth - 1 - i, j, k) = f(hi_idx - pm, j, k);
              f(hi_idx + 1 + i, j, k) = f(bc.depth + pm, j, k);
            } else if constexpr (F == EMFace::Y) {
              const auto pm = j % numInterior;

              f(i, bc.depth - 1 - j, k) = f(i, hi_idx - pm, k);
              f(i, hi_idx + 1 + j, k) = f(i, bc.depth + pm, k);
            } else {
              const auto pm = k % numInterior;

              f(i, j, bc.depth - 1 - k) = f(i, j, hi_idx - pm);
              f(i, j, hi_idx + 1 + k) = f(i, j, bc.depth + pm);
            }
          } // end k-loop
        } // end j-loop
      } // end i-loop
    } // end updateH
  };

  template<EMFace F, EMSide S, Derivative D1, bool Negate>
  struct PMLFunctor {
    using Curl = Diff<D1, Negate>;

    static void apply(auto& f1, const auto& f2, const auto& c1, auto& psi, const auto& b, const auto& c, const size_t i, const size_t j, const size_t k, const size_t x0)
    requires (F == EMFace::X)
    {
      size_t ipml;
      if constexpr (S == EMSide::Lo) { ipml = i; }
      else { ipml = i - x0; }

      psi(ipml, j, k) = b[ipml] * psi(ipml, j, k) + c[ipml] * Curl::apply(f2, i, j, k); //f2(i, j, k) - f2(i - 1, j, k));
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
      else { jpml = j - y0; }

      psi(i, jpml, k) = b[jpml] * psi(i, jpml, k) + c[jpml] * Curl::apply(f2, i, j, k); //(f2(i, j, k) - f2(i, j - 1, k));
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
      else { kpml = k - z0; }

      psi(i, j, kpml) = b[kpml] * psi(i, j, kpml) + c[kpml] * Curl::apply(f2, i, j, k); //(f2(i, j, k) - f2(i, j, k - 1));
      if constexpr (Negate) {
        f1(i, j, k) -= c1(i, j, k) * psi(i, j, kpml);
      } else {
        f1(i, j, k) += c1(i, j, k) * psi(i, j, kpml);
      }
    }
  };

  template<EMFace F, EMSide S, Derivative D1, bool Negate>
  struct BCIntegrator3D {
    using offset_t = tf::electromagnetics::types::IntegratorOffsets;

    static auto apply(auto& f1, const auto& f2, const auto& c1, auto& bc, const offset_t& o) {
      size_t pml_offset;
      if constexpr      (F == EMFace::X) { pml_offset = o.x0; }
      else if constexpr (F == EMFace::Y) { pml_offset = o.y0; }
      else                               { pml_offset = o.z0; }
      if constexpr (Negate) { pml_offset -= 1; } // Correction for H-fields. todo: this will need to be changed when switching E-H stuff

#pragma omp parallel for collapse(3) num_threads(NTHREADS_BC) default(shared)
      for (size_t i = o.x0; i < o.x1; ++i) {
        for (size_t j = o.y0; j < o.y1; ++j) {
          for (size_t k = o.z0; k < o.z1; ++k) {
            PMLFunctor<F, S, D1, Negate>::apply(f1, f2, c1, bc.psi, bc.b, bc.c, i, j, k, pml_offset);
          } // end k-loop
        } // end j-loop
      } // end i-loop
    }
  };

  template<EMFace F, EMSide S, bool isH>
  struct PmlOffsets {
    static types::IntegratorOffsets getPmlOffsets(const auto& bc)
    requires (F == EMFace::X)
    {
      // todo: This whole thing can be moved up into the BCData level so its only executed once at initialization.
      const size_t x0 = bc.offsets.x0 - isH + (S == EMSide::Lo);
      const size_t x1 = bc.offsets.x1 + isH - (S != EMSide::Lo);

      return {x0, x1, bc.offsets.y0, bc.offsets.y1, bc.offsets.z0, bc.offsets.z1};
    }

    static types::IntegratorOffsets getPmlOffsets(const auto& bc)
    requires (F == EMFace::Y)
    {
      // todo: these need to be checked for correctness
      const size_t y0 = bc.offsets.y0 - isH + (S == EMSide::Lo);
      const size_t y1 = bc.offsets.y1 + isH - (S == EMSide::Lo);
      return {bc.offsets.x0, bc.offsets.x1, y0, y1, bc.offsets.z0, bc.offsets.z1};
    }

    static types::IntegratorOffsets getPmlOffsets(const auto& bc)
    requires (F == EMFace::Z)
    {
      // todo: these need to be checked for correctness
      const size_t z0 = bc.offsets.z0 - isH + (S == EMSide::Lo);
      const size_t z1 = bc.offsets.z1 + isH - (S == EMSide::Lo);
      return {bc.offsets.x0, bc.offsets.x1, bc.offsets.y0, bc.offsets.y1, z0, z1};
    }
  };

  template<EMFace F, EMSide S, bool Negate>
  struct Pml3DUpdate : pml_t<F, S> {
    static void updateE(auto& bc, auto& f1, const auto& f2, const auto& c1)
    requires (F == EMFace::X)
    {
      BCIntegrator3D<F, S, Derivative::DX, Negate>::apply(f1, f2, c1, bc, {x0, x1, y0, y1, z0, z1});
    } // end updateE (X-Face)

    static void updateH(auto& bc, auto& f1, const auto& f2, const auto& c1)
    requires (F == EMFace::X)
    {
      BCIntegrator3D<F, S, Derivative::DX, Negate>::apply(f1, f2, c1, bc, {x0, x1, y0, y1, z0, z1});
    } // end updateH (X-Face)

    static void updateE(auto& bc, auto& f1, const auto& f2, const auto& c1)
    requires (F == EMFace::Y)
    {
            PMLFunctor<F, S, Derivative::DY, Negate>::apply(f1, f2, c1, bc.psi, bc.b, bc.c, i, j, k, y0);
    } // end updateE (Y-Face)

    static void updateH(auto& bc, auto& f1, const auto& f2, const auto& c1)
    requires (F == EMFace::Y)
    {
            PMLFunctor<F, S, Derivative::DY, Negate>::apply(f1, f2, c1, bc.psi, bc.b, bc.c, i, j, k, y0 - 1);
    } // end updateH (Y-Face)

    static void updateE(auto& bc, auto& f1, const auto& f2, const auto& c1)
    requires (F == EMFace::Z)
    {
            PMLFunctor<F, S, Derivative::DZ, Negate>::apply(f1, f2, c1, bc.psi, bc.b, bc.c, i, j, k, z0);
    } // end updateE (Z-Face)

    static void updateH(auto& bc, auto& f1, const auto& f2, const auto& c1)
    requires (F == EMFace::Z)
    {
            PMLFunctor<F, S, Derivative::DZ, Negate>::apply(f1, f2, c1, bc.psi, bc.b, bc.c, i, j, k, z0 - 1);
    } // end updateH (Z-Face)
  };
} // end namespace tf::electromagnetics::boundaries

#endif //BOUNDARIES_H
