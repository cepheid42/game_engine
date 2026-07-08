#ifndef CURRENT_DEPOSITION_HPP
#define CURRENT_DEPOSITION_HPP

#include "interpolation.hpp"
#include "particles.hpp"
#include "vec3.hpp"

#include <cassert>
#include <cmath>
#include <print>

namespace tf::particles {
struct CurrentDeposition {
   template<typename Strategy>
   static void deposit(auto& J,
                       const auto p0, const auto p1,
                       const auto& shapeI0, const auto& shapeJ0,
                       const auto& shapeDI, const auto& shapeDJ, const auto& shapeDK,
                       const auto ci, const auto cj, const auto ck, const auto qA)
   {
      using Outer = Strategy::OuterShape;
      using Middle = Strategy::MiddleShape;
      using Inner = Strategy::InnerShape;
      // Return if particle has not moved in this direction (therefore no current)
      if (p0 == p1) { return; }
      for (auto i = Outer::Begin; i <= Outer::End; ++i) {
         const auto& s0i = shapeI0[i];
         const auto& dsi = shapeDI[i];
         for (auto j = Middle::Begin; j <= Middle::End; ++j) {
            const auto& s0j = shapeJ0[j];
            const auto& dsj = shapeDJ[j];
            const auto tmp = -qA * (s0i * s0j + 0.5 * (dsi * s0j + s0i * dsj) + (1.0 / 3.0) * dsj * dsi);
            auto acc = 0.0;
            for (auto k = Inner::Begin; k <= Inner::End - 1; ++k) {
               acc += shapeDK[k] * tmp;
               const auto& [x, y, z] = Strategy::permute(ci + i, cj + j, ck + k);
               #pragma omp atomic update
               J(x, y, z) += acc;
            } // end for(k)
         } // end for(j)
      } // end for(i)
   } // end deposit()

   static void updateJ(const auto& p, auto& emdata, const auto charge) {
      // Aliases for X/Y/Z shape functions
      using XShape = interp::InterpolationShape<x_collapsed ? 1 : interpolation_order>::Type;
      using YShape = interp::InterpolationShape<y_collapsed ? 1 : interpolation_order>::Type;
      using ZShape = interp::InterpolationShape<z_collapsed ? 1 : interpolation_order>::Type;

      // Interpolation Strategies for each J
      using JxStrategy = interp::InterpolationStrategy<0, YShape, ZShape, XShape>;
      using JyStrategy = interp::InterpolationStrategy<1, ZShape, XShape, YShape>;
      using JzStrategy = interp::InterpolationStrategy<2, XShape, YShape, ZShape>;

      static constexpr auto xf = y_collapsed or z_collapsed ? 2.0 : 1.0;
      static constexpr auto yf = x_collapsed or z_collapsed ? 2.0 : 1.0;
      static constexpr auto zf = x_collapsed or y_collapsed ? 2.0 : 1.0;
      // Precompute constants
      static constexpr auto dtAxy = xf / (dt * dx * dy);
      static constexpr auto dtAxz = yf / (dt * dx * dz);
      static constexpr auto dtAyz = zf / (dt * dy * dz);

      // Offsets for Even/Odd order interpolation
      static constexpr vec3 offset{
         interpolation_order == 2 ? 0.5 : 0.0,
         interpolation_order == 2 ? 0.5 : 0.0,
         interpolation_order == 2 ? 0.5 : 0.0
      };

      // Early return if Jdep isn't needed
      if (p.is_disabled()) { return; }

      // Current Density coefficients
      const auto x_coeff = static_cast<double>(p.weight) * charge * dtAyz;
      const auto y_coeff = static_cast<double>(p.weight) * charge * dtAxz;
      const auto z_coeff = static_cast<double>(p.weight) * charge * dtAxy;

      // Find cell indices and determine first relay point
      const vec3 i0 = getCellIndices<double>(p.old_location - offset);
      const vec3 i1 = getCellIndices<double>(p.location - offset);

      const auto same_idx = is_equal(i0, i1);

      const vec3 relay{
         same_idx[0] ? p.location[0] : std::max(i1[0], i0[0]) + offset[0],
         same_idx[1] ? p.location[1] : std::max(i1[1], i0[1]) + offset[1],
         same_idx[2] ? p.location[2] : std::max(i1[2], i0[2]) + offset[2],
      };

      // Calculate normalized locations for first segment
      const vec3 idx0 = i0.to_uint();
      auto p0 = p.old_location - i0;
      auto p1 = relay - i0;

      // Create shape arrays for first segment
      auto s0i = XShape::shape_array(p0[0]);
      auto s0j = YShape::shape_array(p0[1]);
      auto s0k = ZShape::shape_array(p0[2]);
      auto dsi = XShape::ds_array(p1[0], s0i);
      auto dsj = YShape::ds_array(p1[1], s0j);
      auto dsk = ZShape::ds_array(p1[2], s0k);

      // Deposit to J's
      deposit<JxStrategy>(emdata.Jx, p0[0], p1[0], s0j, s0k, dsj, dsk, dsi, idx0[1], idx0[2], idx0[0], x_coeff); // y-z-x
      deposit<JyStrategy>(emdata.Jy, p0[1], p1[1], s0k, s0i, dsk, dsi, dsj, idx0[2], idx0[0], idx0[1], y_coeff); // z-x-y
      deposit<JzStrategy>(emdata.Jz, p0[2], p1[2], s0i, s0j, dsi, dsj, dsk, idx0[0], idx0[1], idx0[2], z_coeff); // x-y-z

      // If particle has changed assignment cell in any direction, do second deposition step
      if (!(same_idx[0] and same_idx[1] and same_idx[2])) {
         // Calculate normalized locations for second segment
         const vec3 idx1 = i1.to_uint();
         p0 = relay - i1;
         p1 = p.location - i1;

         // Create shape arrays for second segment
         s0i = XShape::shape_array(p0[0]);
         s0j = YShape::shape_array(p0[1]);
         s0k = ZShape::shape_array(p0[2]);
         dsi = XShape::ds_array(p1[0], s0i);
         dsj = YShape::ds_array(p1[1], s0j);
         dsk = ZShape::ds_array(p1[2], s0k);

         // Deposit to J's
         deposit<JxStrategy>(emdata.Jx, p0[0], p1[0], s0j, s0k, dsj, dsk, dsi, idx1[1], idx1[2], idx1[0], x_coeff); // y-z-x
         deposit<JyStrategy>(emdata.Jy, p0[1], p1[1], s0k, s0i, dsk, dsi, dsj, idx1[2], idx1[0], idx1[1], y_coeff); // z-x-y
         deposit<JzStrategy>(emdata.Jz, p0[2], p1[2], s0i, s0j, dsi, dsj, dsk, idx1[0], idx1[1], idx1[2], z_coeff); // x-y-z
      }
   } // end updateJ()
   

   static void advance(auto& g, auto& emdata) requires(jdep_enabled) {
      g.jdep_timer.start_timer();
      if (g.is_photons or g.is_tracer) { return; }

      #pragma omp parallel for num_threads(nThreads)
      for (auto pid = 0zu; pid < g.num_particles(); pid++) {
         updateJ(g.particles[pid], emdata, g.charge);
      }
      g.jdep_timer.stop_timer();
   }

   static void advance(auto&, auto&) requires(!jdep_enabled) {}
}; // end struct CurrentDeposition
} // namespace tf::particles

#endif  // CURRENT_DEPOSITION_HPP
