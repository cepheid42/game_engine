#ifndef CURRENT_DEPOSITION_HPP
#define CURRENT_DEPOSITION_HPP

#include "em_data.hpp"
#include "interpolation.hpp"
#include "particles.hpp"

#include "dbg.h"

#include <cmath>

namespace tf::particles {
struct CurrentDeposition {
   template<int D, typename Strategy>
   static void deposit(auto& J, const auto& p0, const auto& p1, const auto& cids, const auto qA) {
      using OShape = typename Strategy::OuterShape::type;
      using MShape = typename Strategy::MiddleShape::type;
      using IShape = typename Strategy::InnerShape::type;

      static constexpr auto iBegin = OShape::Begin;
      static constexpr auto jBegin = MShape::Begin;
      static constexpr auto kBegin = IShape::Begin;
      static constexpr auto   iEnd = OShape::End;
      static constexpr auto   jEnd = MShape::End;
      static constexpr auto   kEnd = IShape::End;

      // if (p0[D] == p1[D]) { return; }

      const auto& [ci, cj, ck] = interp::rotateOrigin<D>(cids);
      const auto& [x0, y0, z0] = interp::rotateOrigin<D>(p0);
      const auto& [x1, y1, z1] = interp::rotateOrigin<D>(p1);

      const auto shapeI0 = OShape::shape_array(x0);
      const auto shapeJ0 = MShape::shape_array(y0);
      const auto shapeK0 = IShape::shape_array(z0);
      const auto shapeI1 = OShape::shape_array(x1);
      const auto shapeJ1 = MShape::shape_array(y1);
      const auto shapeK1 = IShape::shape_array(z1);

      for (int i = iBegin; i <= iEnd; ++i) {
         const auto& s0i = shapeI0[i - iBegin];
         const auto& s1i = shapeI1[i - iBegin];
         const auto dsi = s1i - s0i;
         for (int j = jBegin; j <= jEnd; ++j) {
            const auto& s0j = shapeJ0[j - jBegin];
            const auto& s1j = shapeJ1[j - jBegin];
            const auto dsj = s1j - s0j;
            const auto tmp = -qA * (s0i * s0j + 0.5_fp * (dsi * s0j + s0i * dsj) + (1.0_fp / 3.0_fp) * dsi * dsj);
            auto acc = 0.0_fp;
            for (int k = kBegin; k <= kEnd - 1; ++k) {
               if constexpr (D == 1) {
                  acc += tmp;
               } else {
                  const auto& s0k = shapeK0[k - kBegin];
                  const auto& s1k = shapeK1[k - kBegin];
                  acc += (s1k - s0k) * tmp;
               }

               const auto [x, y, z] = interp::rotateOrigin<D == 2 ? D : !D>(ci + i, cj + j, ck + k);
               #pragma omp atomic update
               J(x, y, z) += acc;
            } // end for(k)
         } // end for(j)
      } // end for(i)
   } // end deposit()


   template<int Support>
   static vec3<compute_t> findRelayPoint(const auto& i0, const auto& i1, const auto& p1) {
      if constexpr (Support % 2 == 0) {
         return {i0[0] == i1[0] ? p1[0] : std::fmax(i0[0], i1[0]) - 0.5,
                 i0[1] == i1[1] ? p1[1] : std::fmax(i0[1], i1[1]) - 0.5,
                 i0[2] == i1[2] ? p1[2] : std::fmax(i0[2], i1[2]) - 0.5};
      }
      else {
         return {i0[0] == i1[0] ? p1[0] : std::fmax(i0[0], i1[0]),
                 i0[1] == i1[1] ? p1[1] : std::fmax(i0[1], i1[1]),
                 i0[2] == i1[2] ? p1[2] : std::fmax(i0[2], i1[2])};

      }
   }

   static void updateJ(const auto& p, auto& emdata, const auto charge) {
      using AssignmentShape = interp::InterpolationShape<interpolation_order>;
      using NGP = interp::InterpolationShape<0>;
      using JxStrategy = interp::InterpolationStrategy<NGP, AssignmentShape, AssignmentShape>;
      using JyStrategy = interp::InterpolationStrategy<AssignmentShape, AssignmentShape, NGP>;
      using JzStrategy = interp::InterpolationStrategy<AssignmentShape, NGP, AssignmentShape>;


      static constexpr auto support = AssignmentShape::Support;
      static constexpr auto offset = AssignmentShape::Order % 2 == 0 ? 0.5_fp : 0.0_fp;

      static constexpr auto dtAxy = 1.0_fp / (dt * dx * dy);
      // static constexpr auto dtAxz = 1.0_fp / (dt * dx * dz);
      static constexpr auto dtAyz = 1.0_fp / (dt * dy * dz);
      if (p.disabled) { return; }

      const auto x_coeff = p.weight * charge * dtAyz;
      const auto y_coeff = -1.0 * p.weight * charge * p.velocity[1] / (dx * dy * dz);
      const auto z_coeff = p.weight * charge * dtAxy;

      const vec3<std::size_t> i0 = getCellIndices(p.old_location + offset);
      const vec3<std::size_t> i1 = getCellIndices(p.location + offset);
      const auto relay = findRelayPoint<support>(i0, i1, p.location);

      auto p0 = p.old_location - i0.as_type<compute_t>();
      auto p1 = relay - i0.as_type<compute_t>();

      deposit<0, JxStrategy>(emdata.Jx, p0, p1, i0, x_coeff);
      deposit<1, JyStrategy>(emdata.Jy, p0, p1, i0, y_coeff);
      deposit<2, JzStrategy>(emdata.Jz, p0, p1, i0, z_coeff);

      if (i0 != i1) {
         p0 = relay - i1.as_type<compute_t>();
         p1 = p.location - i1.as_type<compute_t>();

         deposit<0, JxStrategy>(emdata.Jx, p0, p1, i1, x_coeff);
         // deposit<1, JStrategy>(emdata.Jy, p0, p1, i1, y_coeff);
         deposit<2, JzStrategy>(emdata.Jz, p0, p1, i1, z_coeff);
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
