#ifndef CURRENT_DEPOSITION_HPP
#define CURRENT_DEPOSITION_HPP

#include "particles.hpp"
#include "em_data.hpp"
#include "interpolation.hpp"
#include "constants.hpp"

#include <cmath>

namespace tf::particles {
   namespace detail {
      template<typename StartFunc, typename EndFunc>
      struct TrajectoryAssignment {
         constexpr explicit TrajectoryAssignment(StartFunc&& start_func, EndFunc&& end_func)
         : startShape(std::move(start_func)), endShape(std::move(end_func))
         {}

         [[nodiscard]] constexpr auto S0(const int grid_point) const {
            return startShape(grid_point);
         }

         [[nodiscard]] constexpr auto S1(const int grid_point) const {
            return endShape(grid_point);
         }

         [[nodiscard]] constexpr auto DS(const int grid_point) const {
            return S1(grid_point) - S0(grid_point);
         }

         const StartFunc startShape;
         const EndFunc endShape;
      };

      void depositJx(auto& J, const auto& p0, const auto& p1, const auto& cids, const auto qA) {
         using JitShape = interp::Jit<interp::TSC>;
         using CachedShape = interp::Jit<interp::TSC>;
         using ReducedShape = interp::Jit<interp::CIC>;
         using JitTrajectory = TrajectoryAssignment<JitShape, JitShape>;
         using ReducedTrajectory = TrajectoryAssignment<ReducedShape, ReducedShape>;
         using CachedTrajectory = TrajectoryAssignment<CachedShape, CachedShape>;
         // (x, y, z) -> (1, 2, 0) -> (y, z, x) | Ex D = 0
         // (x, y, z) -> (2, 0, 1) -> (z, x, y) | Ey D = 1
         // (x, y, z) ->           -> (x, y, z) | Ez D = 2
         static constexpr vec3 b0{-1, 0, -1};
         static constexpr vec3 b1{1, 1, 1};
         static constexpr auto third = constants::third<compute_t>;

         if (p0[0] == p1[0]) { return; } // no deposition in this direction

         dbg("Start depositJx");
         const auto& [ci, cj, ck] = cids;
         const auto& [x0, y0, z0] = p0;
         const auto& [x1, y1, z1] = p1;


         JitTrajectory shapeI(JitShape{x0 - ci}, JitShape{x1 - ci});
         ReducedTrajectory shapeJ(ReducedShape{y0 - cj}, ReducedShape{y1 - cj});
         CachedTrajectory shapeK(CachedShape{z0 - ck}, CachedShape{z1 - ck});

         for (int j = b0[1]; j <= b1[1]; ++j) {
            const auto s0j = shapeJ.S0(j);
            const auto dsj = shapeJ.DS(j);

            for (int k = b0[2]; k <= b1[2]; ++k) {
               const auto s0k= shapeK.S0(k);
               const auto dsk = shapeK.DS(k);

               for (int i = b0[0]; i <= b1[0]; ++i) {
                  const auto dsi = shapeI.DS(i);

                  const auto val = qA * dsi * (s0k * s0j + 0.5_fp * (dsk * s0j + s0k * dsj) + third * dsj * dsk);
                  dbg(ci + i, cj + j, ck + k, val);
                  dbg(shapeI.S0(i), shapeI.S1(i), shapeJ.S0(j), shapeJ.S1(j), shapeK.S0(k), shapeK.S1(k));
                  // dbg(shapeI.S0(i), shapeI.S1(i), val);
                  #pragma omp atomic update
                  J(ci + i, cj + j, ck + k) += val;
               }
            }
         }
         dbg("End depositJx");
      }

      void depositJy(auto& Jy, const auto& p0, const auto& p1, const auto& cids, const auto qA) {
         using CachedShape = interp::Cached<interp::TSC>;
         using CachedTrajectory = TrajectoryAssignment<CachedShape, CachedShape>;
         static constexpr auto third = constants::third<compute_t>;

         if (p0[1] == p1[1]) { return; } // No deposition in this direction

         const CachedTrajectory shapeI(CachedShape{p0[0]}, CachedShape{p1[0]});
         const CachedTrajectory shapeK(CachedShape{p0[2]}, CachedShape{p1[2]});

         const auto& [ci, cj, ck] = cids;
         for (int i = -1; i <= 1; i++) {
            const auto s0i = shapeI.S0(ci + i);
            const auto dsi = shapeI.DS(ci + i);
            for (int k = -1; k <= 1; k++) {
               const auto s0k = shapeK.S0(ck + k);
               const auto dsk = shapeK.DS(ck + k);

               // Does this need a minus sign like the normal version?
               #pragma omp atomic update
               Jy(ci + i, cj, ck + k) += qA * third * (s0i * s0k + 0.5_fp * (dsi * s0k + s0i * dsk) + third * dsk * dsi);
            }
         }
      }

      bool findRelayPoint(auto& p0, auto& p1, const auto& pold, const auto& pnew) {
         bool has_second_trajectory = false;
         for (int d = 0; d < 3; d++) {
            const auto i1 = std::floor(pold[d]);
            const auto i2 = std::floor(pnew[d]);
            p0[d] = pold[d];
            if (i1 == i2) {
               p1[d] = pnew[d];
            } else {
               p1[d] = std::max(i1, i2);
               has_second_trajectory = true;
            }
         }
         return has_second_trajectory;
      }

      void updateJ(const auto& p, auto& emdata, const auto charge) {
         if (p.disabled) { return; }

         const auto x_coeff = p.weight * charge / (Ayz * dt);
         const auto z_coeff = p.weight * charge / (Axy * dt);
         const auto y_coeff = p.weight * charge * p.velocity[1] / cellVolume;

         const auto cids = getCIDs(p.location);

         assertm(cids[0] < Ncx - 1 and cids[0] > 1, "updateJ: Invalid particle cell in x.");
         // assertm(cids[1] < Ncy - 1 and cids[1] > 1, "updateJ: Invalid particle cell in y.");
         assertm(cids[2] < Ncz - 1 and cids[2] > 1, "updateJ: Invalid particle cell in z.");

         auto pold = p.old_location - 0.5_fp;
         auto pnew = p.location - 0.5_fp;

         // because first order in Y
         pold[1] += 0.5_fp;
         pnew[1] += 0.5_fp;

         vec3<compute_t> p0{}, p1{};
         const bool second_trajectory = findRelayPoint(p0, p1, pold, pnew);

         depositJx(emdata.Jx, p0, p1, cids, x_coeff);
         // deposit<2>(emdata.Jz, p0, p1, cids, z_coeff);

         if (second_trajectory) {
            // does p1 need to be changed here?
            depositJx(emdata.Jx, p1, pnew, cids, x_coeff);
            // deposit<2>(emdata.Jz, p1, pnew, cids, z_coeff);
         }

         // depositJy(emdata.Jy, pold, pnew, cids, y_coeff);
      }
   } // end namespace tf::particles::detail

   void deposit_current(auto& emdata, const auto& g) {
      #pragma omp parallel for num_threads(nThreads)
      for (std::size_t pid = 0; pid < g.num_particles(); pid++) {
         detail::updateJ(g.particles[pid], emdata, g.charge);
      }
   }
} // end namepsace tf::particles
#endif //CURRENT_DEPOSITION_HPP
