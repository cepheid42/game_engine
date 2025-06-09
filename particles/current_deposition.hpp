#ifndef CURRENT_DEPOSITION_HPP
#define CURRENT_DEPOSITION_HPP

#include "constants.hpp"
#include "em_data.hpp"
#include "interpolation.hpp"
#include "particles.hpp"

#include "dbg.h"

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
         using CachedTSC = interp::Jit<interp::TSC>;
         using CachedCIC = interp::Jit<interp::CIC>;
         using TSCTrajectory = TrajectoryShapeFunc<CachedTSC, CachedTSC>;
         using CICTrajectory = TrajectoryShapeFunc<CachedCIC, CachedCIC>;

         if constexpr (C == 0 or C == 3) {
            return CICTrajectory(CachedCIC{x0}, CachedCIC{x1});
         } else {
            return TSCTrajectory(CachedTSC{x0}, CachedTSC{x1});
         }
      }
   } // end namespace tf::particles::detail

struct CurrentDeposition {
   template<int D>
   static void deposit(auto& J, const auto& p0, const auto& p1, const auto& cids, const auto qA) {
      static constexpr auto third = constants::third<compute_t>;
      static constexpr vec3 b0 = interp::rotateOrigin<D>(-1, 0, -1);
      static constexpr vec3 b1 = interp::rotateOrigin<D>(1, 0, 1);

      if (p0[D] == p1[D]) {
         return;
      }

      const auto& [ci, cj, ck] = interp::rotateOrigin<D>(cids);
      const auto& [x0, y0, z0] = interp::rotateOrigin<D>(p0);
      const auto& [x1, y1, z1] = interp::rotateOrigin<D>(p1);

      const auto shapeI = detail::makeTrajectoryFunction<D>(x0, x1);
      const auto shapeJ = detail::makeTrajectoryFunction<D + 1>(y0, y1);
      const auto shapeK = detail::makeTrajectoryFunction<D + 2>(z0, z1);

      for (int i = b0[0]; i <= b1[0]; ++i) {
         const auto s0i = shapeI.S0(i);
         const auto dsi = shapeI.S1(i) - s0i;

         for (int j = b0[1]; j <= b1[1]; ++j) {
            const auto s0j = shapeJ.S0(j);
            const auto dsj = shapeJ.S1(j) - s0j;

            const auto tmp = -qA * (s0i * s0j + 0.5_fp * (dsi * s0j + s0i * dsj) + third * dsj * dsi);
            auto acc = 0.0;
            for (int k = b0[2]; k <= b1[2]; ++k) {
               acc += shapeK.DS(k) * tmp;
               const auto [x, y, z] = interp::rotateOrigin<D == 2 ? D : !D>(ci + i, cj + j, ck + k);
               // #pragma omp atomic update
               J(x, y, z) += acc;
            } // end for(k)
         } // end for(j)
      } // end for(i)
   } // end deposit()

   static vec3<compute_t> findRelayPoint(const auto& i0, const auto& i1, const auto& x1) {
      // y uses odd-order assignment function (CIC)
      return {
         i0[0] == i1[0] ? x1[0] : std::max(i0[0], i1[0]),
         i0[1] == i1[1] ? x1[1] : std::max(i0[1], i1[1]),
         i0[2] == i1[2] ? x1[2] : std::max(i0[2], i1[2]),
      };
   }

   static void updateJ(const auto& p, auto& emdata, const auto charge) {
      static constexpr auto dtAxy = 1.0_fp / (dt * dx * dy);
      static constexpr auto dtAxz = 1.0_fp / (dt * dx * dz);
      static constexpr auto dtAyz = 1.0_fp / (dt * dy * dz);
      if (p.disabled) { return; }

      const auto x_coeff = p.weight * charge * dtAyz;
      const auto z_coeff = p.weight * charge * dtAxy;
      const auto y_coeff = p.weight * charge * dtAxz;

      auto old_loc_half = p.old_location + 0.5_fp;
      auto new_loc_half = p.location + 0.5_fp;

      const auto old_hcids = getCellIndices<compute_t>(old_loc_half);
      const auto new_hcids = getCellIndices<compute_t>(new_loc_half);
      const auto relay = findRelayPoint(old_hcids, new_hcids, new_loc_half);

      auto p0 = old_loc_half - old_hcids;
      auto p1 = relay - old_hcids;

      p0[1] -= 0.5_fp;
      p1[1] -= 0.5_fp;

      deposit<0>(emdata.Jx, p0, p1, old_hcids.template as_type<std::size_t>(), x_coeff);
      deposit<1>(emdata.Jy, p0, p1, old_hcids.template as_type<std::size_t>(), y_coeff);
      deposit<2>(emdata.Jz, p0, p1, old_hcids.template as_type<std::size_t>(), z_coeff);

      if (old_hcids != new_hcids) {
         p0 = relay - new_hcids;
         p1 = new_loc_half - new_hcids;

         deposit<0>(emdata.Jx, p0, p1, new_hcids.template as_type<std::size_t>(), x_coeff);
         // deposit<1>(emdata.Jy, p0, p1, new_hcids, y_coeff);
         deposit<2>(emdata.Jz, p0, p1, new_hcids.template as_type<std::size_t>(), z_coeff);
      }
   } // end updateJ()

   static void operator()(auto& emdata, const auto& g) {
      #pragma omp parallel for num_threads(nThreads)
      for (std::size_t pid = 0; pid < g.num_particles(); pid++) {
         updateJ(g.particles[pid], emdata, g.charge);
      }
   }
};
} // namespace tf::particles

#endif  // CURRENT_DEPOSITION_HPP
