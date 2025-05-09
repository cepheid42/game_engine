#ifndef CURRENT_DEPOSITION_HPP
#define CURRENT_DEPOSITION_HPP

#include "math_utils.hpp"
#include "particles.hpp"
#include "em_data.hpp"
#include "bc_data.hpp"
#include "morton.hpp"

#include <array>
#include <cmath>
#include <cassert>

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
      // bool active{false};
   };

   static std::array<Segment, 2> split_trajectory(const Particle& p, const std::array<std::size_t, 3>& cidx1) {
      const auto xold = std::floor(p.old_location[0] / dx + 0.5);
      const auto yold = std::floor(p.old_location[1] / dy + 0.5);
      const auto zold = std::floor(p.old_location[2] / dz + 0.5);
      const auto xnew = std::floor(p.location[0] / dx + 0.5);
      const auto ynew = std::floor(p.location[1] / dy + 0.5);
      const auto znew = std::floor(p.location[2] / dz + 0.5);

      const auto xr = x_range[0] + dx * (std::max(xold, xnew) - 0.5);
      const auto yr = y_range[0] + dy * (std::max(yold, ynew) - 0.5);
      const auto zr = z_range[0] + dz * (std::max(zold, znew) - 0.5);

      const auto ix = static_cast<std::size_t>((p.old_location[0] - x_range[0]) / dx);
      const auto iz = static_cast<std::size_t>((p.old_location[2] - z_range[0]) / dz);

      return {
         Segment{{ix, cidx1[1], iz}, p.old_location, {xr, yr, zr}},
         Segment{cidx1, {xr, yr, zr}, p.location}
      };
   }

   struct CurrentDeposition {
      using emdata_t = electromagnetics::EMData;
      using group_t = ParticleGroup;
      using array_t = std::array<std::size_t, 3>;
      using EMFace = electromagnetics::EMFace;
      using EMSide = electromagnetics::EMSide;

      template<int D>
      static void updateJ(auto& J, const auto& qA,
                          const std::array<double, 2>& ws,
                          const std::array<double, 3>& bs0,
                          const std::array<double, 3>& bs1,
                          const std::array<double, 3>& cs0,
                          const std::array<double, 3>& cs1,
                          const std::array<std::size_t, 3>& idxs,
                          const std::array<int, 6>& bounds)
      {
         static constexpr auto third = 1.0_fp / 3.0_fp;
         static constexpr auto sixth = 1.0_fp / 6.0_fp;

         const auto& [i, j, k] = idxs;
         const auto& [x0, x1, y0, y1, z0, z1] = bounds;

         //       x     y     z         p     q     r
         // Jx: {0, 2, 1, 3, 0, 3} -> {0, 2, 1, 3, 0, 3} (x, y, z)
         // Jy: {0, 3, 0, 1, 0, 3} -> {0, 1, 0, 3, 0, 3} (y, x, z)
         // Jz: {0, 3, 1, 3, 0, 2} -> {0, 2, 0, 3, 1, 3} (z, x, y)
         // Jx/Jz: ws = {as0[0] - as1[0, as1[2] - as0[2]]}
         // Jy:    ws = {vy / c}
         std::size_t x, y, z;
         for (int p = x0; p < x1; ++p) {
            for (int q = y0; q < y1; ++q) {
               for (int r = z0; r < z1; ++r) {
                  if constexpr (D == 0) {
                     x = i + p;
                     y = j + q - 1;
                     z = k + r - 1;
                  }
                  else if constexpr (D == 1) {
                     x = i + q - 1;
                     y = j + p;
                     z = k + r - 1;
                  }
                  else {
                     x = i + q - 1;
                     y = j + r - 1,
                     z = k + p;
                  }

                  const auto wT = third * (bs0[q] * cs0[r] + bs1[q] * cs1[r])
                                + sixth * (bs1[q] * cs0[r] + bs0[q] * cs1[r]);

                  // Jx(i + p, j + q - 1, k + r - 1)
                  // Jy(i + q - 1, j + p, k + r - 1)
                  // Jz(i + q - 1, j + r - 1, k + p)
                  J(x, y, z) += qA * ws[p] * wT;
               } // end for(r)
            } // end for(q)
         } // end for(p)
      } // end updateJ()

      static void update_particle(const Particle& p, emdata_t& emdata, const compute_t charge) {
         const auto [ci, cj, ck] = morton_decode(p.code);
         const auto cx = p.weight * charge / (Ayz * dt);
         const auto cy = p.weight * charge / (Axz * dt);
         const auto cz = p.weight * charge / (Axy * dt);

         for (const auto& [cids, p0, p1]: split_trajectory(p, {ci, cj, ck})) {
            const auto xs0 = quadShapes(p0[0]);
            constexpr std::array ys0 = {1.0, 1.0, 1.0};
            const auto zs0 = quadShapes(p0[2]);

            const auto xs1 = quadShapes(p1[0]);
            constexpr std::array ys1 = {1.0, 1.0, 1.0};
            const auto zs1 = quadShapes(p1[2]);

            updateJ<0>(emdata.Jx, cx, {xs0[0] - xs1[0], xs1[2] - xs0[2]}, ys0, ys1, zs0, zs1, cids, {0, 2, 1, 3, 0, 3});
            updateJ<1>(emdata.Jy, cy, {p.velocity[1] / constants::c<compute_t>, 0.0_fp}, xs0, xs1, zs0, zs1, cids, {0, 1, 0, 3, 0, 3});
            updateJ<2>(emdata.Jz, cz, {zs0[0] - zs1[0], zs1[2] - zs0[2]}, xs0, xs1, ys0, ys1, cids, {0, 2, 0, 3, 1, 3});
         } // end for(trajectory)
      }

      static void operator()(const group_t& g, emdata_t& emdata) {
         // #pragma omp parallel for num_threads(nThreads)
         for (std::size_t pid = 0; pid < g.num_particles(); pid++) {
            update_particle(g.particles[pid], emdata, g.charge);
         }
      }
   }; // end struct CurrentDeposition
} // end namepsace tf::particles
#endif //CURRENT_DEPOSITION_HPP
