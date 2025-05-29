#ifndef CURRENT_DEPOSITION_HPP
#define CURRENT_DEPOSITION_HPP

#include "constants.hpp"
#include "em_data.hpp"
#include "interpolation.hpp"
#include "particles.hpp"

#include <cmath>

namespace tf::particles {
   namespace detail {
      template <typename StartFunc, typename EndFunc>
      struct TrajectoryShapeFunc {
         constexpr explicit TrajectoryShapeFunc(StartFunc&& start_func, EndFunc&& end_func)
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

      template<int C>
      constexpr auto makeTrajectoryFunction(const auto x0, const auto x1) {
         using CachedTSC = interp::Cached<interp::TSC>;
         using CachedCIC = interp::Cached<interp::CIC>;
         using CachedTrajectory = TrajectoryShapeFunc<CachedTSC, CachedTSC>;
         using ReducedTrajectory = TrajectoryShapeFunc<CachedCIC, CachedCIC>;

         if constexpr (C == 0 or C == 3) {
            return ReducedTrajectory(CachedCIC{x0}, CachedCIC{x1});
         } else {
            return CachedTrajectory(CachedTSC{x0}, CachedTSC{x1});
         }
      }

      template<int D>
      void deposit(auto& J, const auto& p0, const auto& p1, const auto& cids, const auto qA) {
         static constexpr auto third = constants::third<compute_t>;
         static constexpr vec3 b0{-1, -1, -1};
         static constexpr vec3 b1 = interp::rotateOrigin<D>(1, D == 1 ? -1 : 0, 1);

         if (p0[D] == p1[D]) { return; }

         const auto& [ci, cj, ck] = interp::rotateOrigin<D>(cids);
         const auto& [x0, y0, z0] = interp::rotateOrigin<D>(p0);
         const auto& [x1, y1, z1] = interp::rotateOrigin<D>(p1);

         const auto shapeI = makeTrajectoryFunction<D>(x0, x1);
         const auto shapeJ = makeTrajectoryFunction<D + 1>(y0, y1);
         const auto shapeK = makeTrajectoryFunction<D + 2>(z0, z1);

         for (int i = b0[0]; i <= b1[0]; ++i) {
            const auto s0i = shapeI.S0(i);
            const auto dsi = shapeI.S1(i) - s0i;

            for (int j = b0[1]; j <= b1[1]; ++j) {
               const auto s0j = shapeJ.S0(j);
               const auto dsj = shapeJ.S1(j) - s0j;

               const auto tmp = qA * (s0i * s0j + 0.5_fp * (dsi * s0j + s0i * dsj) + third * dsj * dsi);
               auto acc = 0.0;
               for (int k = b0[2]; k <= b1[2]; ++k) {
                  acc += shapeK.DS(k) * tmp;

                  const auto [x, y, z] = interp::rotateOrigin<D == 2 ? D : !D>(ci + i, cj + j, ck + k);
                  #pragma omp atomic update
                  J(x, y, z) += acc;
               } // end for(k)
            } // end for(j)
         } // end for(i)
      } // end deposit()

      vec3<compute_t> findRelayPoint(const auto& i0, const auto& i1, const auto& x1) {
         // y uses even-order assignment function (CIC)
         return {
            i0[0] == i1[0] ? x1[0] : 0.5_fp * (i0[0] + i1[0]),
            i0[1] == i1[1] ? x1[1] : std::max(i0[1], i1[1]),
            i0[2] == i1[2] ? x1[2] : 0.5_fp * (i0[2] + i1[2]),
         };
      }

      vec3<compute_t> getOffsets(const auto& loc) {
         // y-offset is not shifted by 0.5
         return {std::floor(loc[0] + 0.5_fp),
                 std::floor(loc[1]),
                 std::floor(loc[2] + 0.5_fp)};
      }

      void updateJ(const auto& p, auto& emdata, const auto charge) {
         if (p.disabled) { return; }

         const auto x_coeff = p.weight * charge * dtAyz;
         const auto z_coeff = p.weight * charge * dtAxy;
         const auto y_coeff = p.weight * charge * dtAxz;

         auto cids = getCIDs(p.old_location + 0.5_fp);
         const auto i0 = getOffsets(p.old_location);
         const auto i1 = getOffsets(p.location);
         const auto relay = findRelayPoint(i0, i1, p.location);

         auto p0 = p.old_location - i0;
         auto p1 = relay - i0;

         deposit<0>(emdata.Jx, p0, p1, cids, x_coeff);
         deposit<1>(emdata.Jy, p0, p1, cids, y_coeff);
         deposit<2>(emdata.Jz, p0, p1, cids, z_coeff);

         if (i0 != i1) {
            p0 = relay - i1;
            p1 = p.location - i1;
            for (int d = 0; d < 3; ++d) {
               cids[d] += static_cast<int>(i1[d] - i0[d]);
            }

            deposit<0>(emdata.Jx, p0, p1, cids, x_coeff);
            // deposit<1>(emdata.Jy, p0, p1, cids, y_coeff);
            deposit<2>(emdata.Jz, p0, p1, cids, z_coeff);
         }
      } // end updateJ()
   } // namespace detail

   void deposit_current(auto& emdata, const auto& g) {
      #pragma omp parallel for num_threads(nThreads)
      for (std::size_t pid = 0; pid < g.num_particles(); pid++) {
         detail::updateJ(g.particles[pid], emdata, g.charge);
      }
   }
} // namespace tf::particles

#endif  // CURRENT_DEPOSITION_HPP
