#ifndef CURRENT_DEPOSITION_HPP
#define CURRENT_DEPOSITION_HPP

#include "interpolation.hpp"
#include "particles.hpp"
#include "vec3.hpp"

#include <cmath>

namespace tf::particles {
struct CurrentDeposition {
   template<int D, typename Shape>
   static void deposit(auto& J,
                       const auto p0, const auto p1,
                       const auto& shapeI0, const auto& shapeJ0,
                       const auto& shapeDI, const auto& shapeDJ, const auto& shapeDK,
                       const auto ci, const auto cj, const auto ck, const auto qA)
   {
      if (p0 == p1) { return; }

      for (int i = Shape::Begin; i <= Shape::End; ++i) {
         const auto& s0i = shapeI0[i - Shape::Begin];
         const auto& dsi = shapeDI[i - Shape::Begin];
         for (int j = Shape::Begin; j <= Shape::End; ++j) {
            const auto& s0j = shapeJ0[j - Shape::Begin];
            const auto& dsj = shapeDJ[j - Shape::Begin];
            const auto tmp = -qA * (s0i * s0j + 0.5 * (dsi * s0j + s0i * dsj) + (1.0 / 3.0) * dsj * dsi);
            auto acc = 0.0;
            for (int k = Shape::Begin; k <= Shape::End - 1; ++k) {
               acc += shapeDK[k - Shape::Begin] * tmp;
               const auto [x, y, z] = interp::rotateOrigin<D == 2 ? D : !D>(ci + i, cj + j, ck + k);
               #pragma omp atomic update
               J(x, y, z) += acc;
            } // end for(k)
         } // end for(j)
      } // end for(i)
   } // end deposit()


   template<int D, typename Shape>
   static void deposit2d(auto& J,
                         const auto p0, const auto p1,
                         const auto& shapeI0, const auto& shapeJ0,
                         const auto& shapeDI, const auto& shapeDJ, const auto&,
                         const auto ci, const auto cj, const auto, const auto qA)
   {
      if (p0 == p1) { return; }

      for (int i = Shape::Begin; i <= Shape::End; ++i) {
         const auto& s0i = shapeI0[i - Shape::Begin];
         const auto& dsi = shapeDI[i - Shape::Begin];
         for (int j = Shape::Begin; j <= Shape::End; ++j) {
            const auto& s0j = shapeJ0[j - Shape::Begin];
            const auto& dsj = shapeDJ[j - Shape::Begin];
            const auto result = qA * (s0i * s0j + 0.5 * (dsi * s0j + s0i * dsj) + (1.0 / 3.0) * dsj * dsi);
            const auto [x, y, z] = interp::rotateOrigin<D == 2 ? D : !D>(ci + i, cj + j, 0lu);
            #pragma omp atomic update
            J(x, y, z) += result;
         } // end for(j)
      } // end for(i)
   } // end deposit()


   static void updateJ(const auto& p, auto& emdata, const auto charge) {
      using Shape = interp::InterpolationShape<interpolation_order>::Type;

      static constexpr auto dtAxy = 1.0 / (dt * dx * dy);
      // static constexpr auto dtAxz = 1.0 / (dt * dx * dz);
      static constexpr auto dtAyz = 1.0 / (dt * dy * dz);

      if (p.disabled) { return; }

      const auto x_coeff = p.weight * charge * dtAyz;
      // const auto y_coeff = p.weight * charge * dtAxz;
      const auto y_coeff = p.weight * charge * p.velocity[1] / (dx * dy * dz);
      const auto z_coeff = p.weight * charge * dtAxy;

      static constexpr auto offset = Shape::Order % 2 == 0 ? 0.5 : 1.0;
      const vec3 i0 = getCellIndices<double>(p.old_location + offset);
      const vec3 i1 = getCellIndices<double>(p.location + offset);

      const vec3 idx0 = i0.template as_type<std::size_t>();
      const vec3 idx1 = i1.template as_type<std::size_t>();

      const vec3<double> relay{
         i1[0] == i0[0] ? p.location[0] : std::fmax(i1[0], i0[0]) - offset,
         i1[1] == i0[1] ? p.location[1] : std::fmax(i1[1], i0[1]) - offset,
         i1[2] == i0[2] ? p.location[2] : std::fmax(i1[2], i0[2]) - offset,
      };

      auto p0 = p.old_location - i0;
      auto p1 = relay - i0;

      auto s0i = Shape::shape_array(p0[0]);
      auto s0j = Shape::shape_array(p0[1]);
      auto s0k = Shape::shape_array(p0[2]);
      auto dsi = Shape::ds_array(p1[0], s0i);
      auto dsj = Shape::ds_array(p1[1], s0j);
      auto dsk = Shape::ds_array(p1[2], s0k);

      deposit<0, Shape>(emdata.Jx, p0[0], p1[0], s0j, s0k, dsj, dsk, dsi, idx0[1], idx0[2], idx0[0], x_coeff); // y-z-x
      // deposit<1, Shape>(emdata.Jy, p0[1], p1[1], s0k, s0i, dsk, dsi, dsj, idx0[2], idx0[0], idx0[1], y_coeff); // z-x-y
      deposit2d<1, Shape>(emdata.Jy, p0[1], p1[1], s0k, s0i, dsk, dsi, dsj, idx0[2], idx0[0], idx0[1], y_coeff); // z-x-y
      deposit<2, Shape>(emdata.Jz, p0[2], p1[2], s0i, s0j, dsi, dsj, dsk, idx0[0], idx0[1], idx0[2], z_coeff); // x-y-z

      if (i0 != i1) {
         p0 = relay - i1;
         p1 = p.location - i1;
         s0i = Shape::shape_array(p0[0]);
         s0j = Shape::shape_array(p0[1]);
         s0k = Shape::shape_array(p0[2]);
         dsi = Shape::ds_array(p1[0], s0i);
         dsj = Shape::ds_array(p1[1], s0j);
         dsk = Shape::ds_array(p1[2], s0k);

         deposit<0, Shape>(emdata.Jx, p0[0], p1[0], s0j, s0k, dsj, dsk, dsi, idx1[1], idx1[2], idx1[0], x_coeff); // y-z-x
         // deposit<1, Shape>(emdata.Jy, p0[1], p1[1], s0k, s0i, dsk, dsi, dsj, idx1[2], idx1[0], idx1[1], y_coeff); // z-x-y
         deposit<2, Shape>(emdata.Jz, p0[2], p1[2], s0i, s0j, dsi, dsj, dsk, idx1[0], idx1[1], idx1[2], z_coeff); // x-y-z
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
