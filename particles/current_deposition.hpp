#ifndef CURRENT_DEPOSITION_HPP
#define CURRENT_DEPOSITION_HPP

#include "math_utils.hpp"
#include "particles.hpp"
#include "em_data.hpp"
#include "../electromagnetics/bc_data.hpp"
#include "morton.hpp"

#include <array>
#include <cmath>

namespace tf::particles {
   static std::array<compute_t, 3> quadShapes(const compute_t v) {
      return {
         0.5_fp * math::SQR(0.5_fp - v),
         0.75_fp - math::SQR(v),
         0.5_fp * math::SQR(0.5_fp + v)
      };
   }

   struct Segment {
      std::array<std::size_t, 3> cids = {0, 0, 0};
      vec3<compute_t> p0{};
      vec3<compute_t> p1{};
      bool active{false};
   };

   static std::array<Segment, 2> split_trajectory(const Particle& p, const std::array<std::size_t, 3>& cidx1) {
      const auto xoff = std::floor(p.old_location[0]);
      const auto yoff = std::floor(p.old_location[1]);
      const auto zoff = std::floor(p.old_location[2]);

      if (xoff == 0.0 and yoff == 0.0 and zoff == 0.0) {
         // Particle did not leave current cell
         return {
            Segment{cidx1, p.old_location, p.location, true},
            Segment{cidx1, p.old_location, p.location, false}
         };
      }

      const auto xr0 = xoff == 0.0_fp ? p.location[0] : std::max(0.0_fp, -xoff);
      const auto yr0 = yoff == 0.0_fp ? p.location[1] : std::max(0.0_fp, -yoff);
      const auto zr0 = zoff == 0.0_fp ? p.location[2] : std::max(0.0_fp, -zoff);

      const vec3 pr0{xr0, yr0, zr0};
      const vec3 pr1{xr0 + xoff, yr0 + yoff, zr0 + zoff};

      // todo: pray to the gods that none of the cell id's are zero...
      const auto xo = static_cast<int>(xoff);
      const auto yo = static_cast<int>(yoff);
      const auto zo = static_cast<int>(zoff);
      // Have to offset old_location so that p0 is normalized to the old cell, since p1 already is
      // e.g. old={-1.5,...} -> pr0[{1.0,...} instead of pold={0.5,...} -> pr0{1.0,...}
      const vec3 pold{p.old_location[0] - xo, p.old_location[1] - yo, p.old_location[2] - zo};
      const std::array cidx0 = {cidx1[0] + xo, cidx1[1] + yo, cidx1[2] + zo};
      return {
         Segment{cidx0, pold, pr0, true},
         Segment{cidx1, pr1, p.location, true}
      };
   }

   struct CurrentDeposition {
      using emdata_t = electromagnetics::EMData;
      using group_t = ParticleGroup;
      using array_t = std::array<std::size_t, 3>;
      using EMFace = electromagnetics::EMFace;
      using EMSide = electromagnetics::EMSide;

      template<int D>
      static void updateJ(auto& J, const auto& xs0, const auto& xs1, const auto& ys0, const auto& ys1, const auto& zs0,
                          const auto& zs1, const auto& qA, const std::array<std::size_t, 3>& idxs,
                          const std::array<std::size_t, 3>& bounds) {
         static constexpr auto third = 1.0_fp / 3.0_fp;
         static constexpr auto sixth = 1.0_fp / 6.0_fp;

         const auto& [i, j, k] = idxs;
         const auto& [i1, j1, k1] = bounds;

         if constexpr (D == 0) {
            const std::array<compute_t, 2> ws = {xs0[0] - xs1[0], xs1[2] - xs0[2]};
            for (std::size_t p = 0; p < i1; p++) {
               for (std::size_t q = 0; q < j1; q++) {
                  for (std::size_t r = 0; r < k1; r++) {
                     const auto wT = third * (ys0[q] * zs0[r] + ys1[q] * zs1[r])
                                   + sixth * (ys1[q] * zs0[r] + ys0[q] * zs1[r]);
                     #pragma omp atomic update
                     J(i + p, j + q, k + r) += qA * ws[p] * wT;
                  } // end for(p)
               } // end for(q)
            } // end for(r)
         } // end if(D == 0)
         else if constexpr (D == 1) {
            const std::array<compute_t, 2> ws = {ys0[0] - ys1[0], ys1[2] - ys0[2]};
            for (std::size_t p = 0; p < i1; p++) {
               for (std::size_t q = 0; q < j1; q++) {
                  for (std::size_t r = 0; r < k1; r++) {
                     const auto wT = third * (xs0[p] * zs0[r] + xs1[p] * zs1[r])
                                   + sixth * (xs1[p] * zs0[r] + xs0[p] * zs1[r]);
                     #pragma omp atomic update
                     J(i + p, j + q, k + r) += qA * ws[q] * wT;
                  } // end for(p)
               } // end for(q)
            } // end for(r)
         } // end if(D == 1)
         else {
            const std::array<compute_t, 2> ws = {zs0[0] - zs1[0], zs1[2] - zs0[2]};
            for (std::size_t p = 0; p < i1; p++) {
               for (std::size_t q = 0; q < j1; q++) {
                  for (std::size_t r = 0; r < k1; r++) {
                     const auto wT = third * (xs0[p] * ys0[q] + xs1[p] * ys1[q])
                                   + sixth * (xs1[p] * ys0[q] + xs0[q] * ys1[q]);
                     #pragma omp atomic update
                     J(i + p, j + q, k + r) += qA * ws[r] * wT;
                  } // end for(p)
               } // end for(q)
            } // end for(r)
         } // end if
      } // end updateJ()

      static void update_particle(const Particle& p, emdata_t& emdata, const compute_t charge) {
         const auto [ci, cj, ck] = morton_decode(p.code);
         for (const auto& [cids, p0, p1, active]: split_trajectory(p, {ci, cj, ck})) {
            if (!active) { continue; }
            const auto& [i, j, k] = cids;

            const auto xs0 = quadShapes(p0[0]);
            const auto ys0 = quadShapes(p0[1]);
            const auto zs0 = quadShapes(p0[2]);

            const auto xs1 = quadShapes(p1[0]);
            const auto ys1 = quadShapes(p1[1]);
            const auto zs1 = quadShapes(p1[2]);

            const auto cx = p.weight * charge / (Ayz * dt);
            const auto cy = p.weight * charge / (Axz * dt);
            const auto cz = p.weight * charge / (Axy * dt);

            // todo: j is fixed at 0 for 2d, 3d will require j-1 in x and z
            updateJ<0>(emdata.Jx, xs0, xs1, ys0, ys1, zs0, zs1, cx, {    i, j, k - 1}, {2, 1, 3});
            updateJ<1>(emdata.Jy, xs0, xs1, ys0, ys1, zs0, zs1, cy, {i - 1, j, k - 1}, {3, 1, 3});
            updateJ<2>(emdata.Jz, xs0, xs1, ys0, ys1, zs0, zs1, cz, {i - 1, j, k    }, {3, 1, 2});
         } // end for(trajectory)
      }

      static void operator()(const group_t& g, emdata_t& emdata) {
         #pragma omp parallel for num_threads(nThreads)
         for (std::size_t pid = 0; pid < g.num_particles(); pid++) {
            update_particle(g.particles[pid], emdata, g.charge);
         }
      }
   }; // end struct CurrentDeposition
} // end namepsace tf::particles
#endif //CURRENT_DEPOSITION_HPP
