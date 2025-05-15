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

   template<int D>
   auto rotateOrigin(const auto& p) {
      if constexpr (D == 0) {
         return decltype(p){p[1], p[2], p[0]};
      }
      else if constexpr (D == 1) {
         return decltype(p){p[2], p[0], p[1]};
      }
      else {
         return p;
      }
   }

   struct CurrentDeposition {
      using emdata_t = electromagnetics::EMData;
      using group_t = ParticleGroup;
      using array_t = std::array<std::size_t, 3>;
      using EMFace = electromagnetics::EMFace;
      using EMSide = electromagnetics::EMSide;

      static constexpr compute_t third = 1.0_fp / 3.0_fp;

      static auto shape(auto x) {
         const auto absx = std::abs(x);
         if (absx < 0.5) {
            return 0.75 - math::SQR(x);
         }
         return 0.5 * math::SQR(1.5 - absx);
      }

      template<int D>
      static void deposit(auto& J, const auto& p0, const auto& p1, const auto& cids, const auto qA) {
         static constexpr int i0 = D == 0 ? 0 : -1;
         static constexpr int i1 = 1;

         static constexpr int j0 = -1;
         static constexpr int j1 = 1;

         static constexpr int k0 = -1;
         static constexpr int k1 = 1;

         if (p0[D] == p1[D]) { return; /* no deposition in this direction */ }

         // rotate elements to match loop structure
         const auto [ci, cj, ck] = rotateOrigin<D>(cids);
         const auto [x0, y0, z0] = rotateOrigin<D>(p0);
         const auto [x1, y1, z1] = rotateOrigin<D>(p1);

         for (int i = i0; i <= i1; ++i) {
            const auto s0i = shape(x0 - i);
            const auto dsi = shape(x1 - i) - s0i;

            for (int j = j0; j <= j1; ++j) {
               const auto s0j = shape(y0 - j);
               const auto dsj = shape(y1 - j) - s0j;

               for (int k = k0; k <= k1; ++k) {
                  const auto s0k = shape(z0 - k);
                  const auto dsk = shape(z1 - k) - s0k;

                  // undo the rotation to get proper indices back
                  const auto [x, y, z] = rotateOrigin<D == 2 ? D : !D>(vec3{i + ci, j + cj, k + ck});
                  J(x, y, z) += -qA * dsk * (s0i * s0j + 0.5_fp * (dsi * s0j + s0i * dsj) + third * dsj * dsi);
               }
            }
         }
      }

      static void depositJy(auto& Jy, const auto& p0, const auto& p1, const auto& cids, const auto qA) {
         if (p0[1] == p1[1]) { return; }

         const auto& [ci, cj, ck] = cids;
         for (int i = -1; i <= 1; i++) {
            const auto s0i = shape(p0[0] - i);
            const auto dsi = shape(p1[0] - i) - s0i;
            for (int k = -1; k <= 1; k++) {
               const auto s0k = shape(p0[2] - k);
               const auto dsk = shape(p1[2] - k) - s0k;

               Jy(ci + i, cj, ck + k) += qA * third * (s0i * s0k + 0.5_fp * (dsi * s0k + s0i * dsk) + third * dsk * dsi);
            }
         }
      }

      static void updateJ(const Particle& p, emdata_t& emdata, const compute_t charge) {
         const auto cids = morton_decode(p.code);

         const vec3 off = {std::floor(p.old_location[0]),
                           std::floor(p.old_location[1]),
                           std::floor(p.old_location[2])};

         assert(off[1] == 0);

         vec3<compute_t> i1{}, i2{}, p0{}, p1{};

         for (int d = 0; d < 3; d++) {
            i1[d] = std::floor(p.old_location[d] + 0.5_fp);
            i2[d] = std::floor(p.location[d] + 0.5_fp);

            p0[d] = p.old_location[d] - off[d];
            p1[d] = (i1 == i2) ? p.location[0] : 0.5_fp * (i1[d] + i2[d]);
         }

         deposit<0>(emdata.Jx, p0, p1, cids, p.weight * charge / (Ayz * dt));
         deposit<2>(emdata.Jz, p0, p1, cids, p.weight * charge / (Axy * dt));

         if (i1 != i2) {
            deposit<0>(emdata.Jx, p1, p.location, cids, p.weight * charge / (Ayz * dt));
            deposit<2>(emdata.Jz, p1, p.location, cids, p.weight * charge / (Axy * dt));
         }

         depositJy(emdata.Jy, p.old_location - off, p.location, cids, p.weight * charge * p.velocity[1] / cellVolume);
      }

      static void operator()(const group_t& g, emdata_t& emdata) {
         // #pragma omp parallel for num_threads(nThreads)
         for (std::size_t pid = 0; pid < g.num_particles(); pid++) {
            updateJ(g.particles[pid], emdata, g.charge);
         }
      }
   }; // end struct CurrentDeposition
} // end namepsace tf::particles
#endif //CURRENT_DEPOSITION_HPP
