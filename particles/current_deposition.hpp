#ifndef CURRENT_DEPOSITION_HPP
#define CURRENT_DEPOSITION_HPP

#include "interpolation.hpp"
#include "particles.hpp"
#include "vec3.hpp"

#include <cmath>

namespace tf::particles {
struct CurrentDeposition {
   template<int D, typename Strategy>
   static void deposit(auto& J,
                       const auto p0, const auto p1,
                       const auto& shapeI0, const auto& shapeJ0,
                       const auto& shapeDI, const auto& shapeDJ, const auto& shapeDK,
                       const auto ci, const auto cj, const auto ck, const auto qA)
   {
      using Outer = typename Strategy::OuterShape;
      using Middle = typename Strategy::MiddleShape;
      using Inner = typename Strategy::InnerShape;
      // Return if particle has not moved in this direction (therefore no current)
      if (p0 == p1) { return; }
      for (int i = Outer::Begin; i <= Outer::End; ++i) {
         const auto idx = i - Outer::Begin;
         const auto& s0i = shapeI0[idx];
         const auto& dsi = shapeDI[idx];
         for (int j = Middle::Begin; j <= Middle::End; ++j) {
            const auto jdx = j - Middle::Begin;
            const auto& s0j = shapeJ0[jdx];
            const auto& dsj = shapeDJ[jdx];
            const auto tmp = -qA * (s0i * s0j + 0.5 * (dsi * s0j + s0i * dsj) + (1.0 / 3.0) * dsj * dsi);
            // Ask Ayden why this accumulator works if you want to know
            auto acc = 0.0;
            for (int k = Inner::Begin; k <= Inner::End - 1; ++k) {
               acc += shapeDK[k - Inner::Begin] * tmp;
               const auto [x, y, z] = interp::rotateOrigin<D == 2 ? D : !D>(ci + i, cj + j, ck + k);
               #pragma omp atomic update
               J(x, y, z) += acc;
            } // end for(k)
         } // end for(j)
      } // end for(i)
   } // end deposit()


   template<int D, typename Strategy>
   requires(is_2D_XZ and D == 1)
   static void deposit(auto& J,
                       const auto p0, const auto p1,
                       const auto& shapeI0, const auto& shapeJ0,
                       const auto& shapeDI, const auto& shapeDJ, const auto&,
                       const auto ci, const auto cj, const auto, const auto qA)
   {
      // This is the 2D version of the deposition, intended mainly for the Y-direction currently
      // It only interpolates to outer loops (x/z) and the charge density qA includes the velocity
      // and is positive (not -qA like other version). See EZ-PIC or Esirkepov for more info.
      using Outer = typename Strategy::OuterShape;
      using Middle = typename Strategy::MiddleShape;
      // Return if particle has not moved in this direction (therefore no current)
      if (p0 == p1) { return; }
      for (int i = Outer::Begin; i <= Outer::End; ++i) {
         const auto idx = i - Outer::Begin;;
         const auto& s0i = shapeI0[idx];
         const auto& dsi = shapeDI[idx];
         for (int j = Middle::Begin; j <= Middle::End; ++j) {
            const auto jdx = j - Middle::Begin;
            const auto& s0j = shapeJ0[jdx];
            const auto& dsj = shapeDJ[jdx];
            const auto [x, y, z] = interp::rotateOrigin<D == 2 ? D : !D>(ci + i, cj + j, 0lu);
            #pragma omp atomic update
            J(x, y, z) +=  qA * (s0i * s0j + 0.5 * (dsi * s0j + s0i * dsj) + (1.0 / 3.0) * dsj * dsi);
         } // end for(j)
      } // end for(i)
   } // end deposit()


   static void updateJ(const auto& p, auto& emdata, const auto charge) {
      // Aliases for X/Y/Z shape functions
      using XShape = interp::InterpolationShape<interpolation_order>::Type;
      using YShape = interp::InterpolationShape<is_2D_XZ ? 1 : interpolation_order>::Type;
      using ZShape = interp::InterpolationShape<interpolation_order>::Type;
      // Interpolation Strategies for each J
      using JxStrategy = interp::InterpolationStrategy<YShape, ZShape, XShape>;
      using JyStrategy = interp::InterpolationStrategy<ZShape, XShape, YShape>;
      using JzStrategy = interp::InterpolationStrategy<XShape, YShape, ZShape>;
      // Precompute constants
      static constexpr auto dtAxy = 1.0 / (dt * dx * dy);
      static constexpr auto dtAxz = is_2D_XZ ? 1.0 / (dx * dy * dz) : 1.0 / (dt * dx * dz);
      static constexpr auto dtAyz = 1.0 / (dt * dy * dz);
      auto y_vel = 1.0;
      if constexpr (is_2D_XZ) { y_vel = p.velocity[1]; }
      // Early return if Jdep isn't needed
      if (p.disabled) { return; }
      // Current Density coefficients
      const auto x_coeff = p.weight * charge * dtAyz;
      const auto y_coeff = p.weight * charge * dtAxz * y_vel;
      const auto z_coeff = p.weight * charge * dtAxy;
      // Offsets for Even/Odd order interpolation
      static constexpr vec3 offsets{
         XShape::Order % 2 == 0 ? 0.5 : 1.0,
         YShape::Order % 2 == 0 ? 0.5 : 1.0,
         ZShape::Order % 2 == 0 ? 0.5 : 1.0,
      };
      // Find cell indices and determine first relay point
      const vec3<double> i0 = getCellIndices<double>(p.old_location + offsets);
      const vec3<double> i1 = getCellIndices<double>(p.location + offsets);
      const auto same_idx = is_equal(i0, i1);
      const vec3<double> relay{
         same_idx[0] ? p.location[0] : std::fmax(i1[0], i0[0]) - offsets[0],
         same_idx[1] ? p.location[1] : std::fmax(i1[1], i0[1]) - offsets[1],
         same_idx[2] ? p.location[2] : std::fmax(i1[2], i0[2]) - offsets[2],
      };
      // Calculate normalized locations for first segment
      auto p0 = p.old_location - i0;
      auto p1 = relay - i0;
      const vec3 idx0 = i0.as_type<std::size_t>();
      const vec3 idx1 = i1.as_type<std::size_t>();
      // Create shape arrays for first segment
      auto s0i = XShape::shape_array(p0[0]);
      auto s0j = YShape::shape_array(p0[1]);
      auto s0k = ZShape::shape_array(p0[2]);
      auto dsi = XShape::ds_array(p1[0], s0i);
      auto dsj = YShape::ds_array(p1[1], s0j);
      auto dsk = ZShape::ds_array(p1[2], s0k);
      // Deposit to J's
      deposit<0, JxStrategy>(emdata.Jx, p0[0], p1[0], s0j, s0k, dsj, dsk, dsi, idx0[1], idx0[2], idx0[0], x_coeff); // y-z-x
      deposit<1, JyStrategy>(emdata.Jy, p0[1], p1[1], s0k, s0i, dsk, dsi, dsj, idx0[2], idx0[0], idx0[1], y_coeff); // z-x-y
      deposit<2, JzStrategy>(emdata.Jz, p0[2], p1[2], s0i, s0j, dsi, dsj, dsk, idx0[0], idx0[1], idx0[2], z_coeff); // x-y-z
      // If particle has changed assignment cell in any direction, do second deposition step
      if (!(same_idx[0] and same_idx[1] and same_idx[2])) {
         // Calculate normalized locations for second segment
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
         deposit<0, JxStrategy>(emdata.Jx, p0[0], p1[0], s0j, s0k, dsj, dsk, dsi, idx1[1], idx1[2], idx1[0], x_coeff); // y-z-x
         deposit<1, JyStrategy>(emdata.Jy, p0[1], p1[1], s0k, s0i, dsk, dsi, dsj, idx1[2], idx1[0], idx1[1], y_coeff); // z-x-y
         deposit<2, JzStrategy>(emdata.Jz, p0[2], p1[2], s0i, s0j, dsi, dsj, dsk, idx1[0], idx1[1], idx1[2], z_coeff); // x-y-z
      }
   } // end updateJ()
   

   static void advance(const auto& g, auto& emdata) {
      #pragma omp parallel for num_threads(nThreads)
      for (std::size_t pid = 0; pid < g.num_particles(); pid++) {
         updateJ(g.particles[pid], emdata, g.charge);
      }
   }
}; // end struct CurrentDeposition
} // namespace tf::particles

#endif  // CURRENT_DEPOSITION_HPP
