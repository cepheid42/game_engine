#ifndef CURRENT_DEPOSITION_HPP
#define CURRENT_DEPOSITION_HPP

// #include "em_data.hpp"
#include "interpolation.hpp"
#include "particles.hpp"


#include <cmath>

namespace tf::particles {
struct CurrentDeposition {
   template<int D, typename Shape>
   static void deposit(auto& J,
                       const auto& p0, const auto& p1,
                       const auto& shapeI0, const auto& shapeJ0, const auto& shapeK0,
                       const auto& shapeI1, const auto& shapeJ1, const auto& shapeK1,
                       const auto& cids, const auto qA)
   {
      static constexpr auto Begin = Shape::Begin;
      static constexpr auto   End = Shape::End;
      
      if (p0[D] == p1[D]) { return; }

      const auto& [ci, cj, ck] = interp::rotateOrigin<D>(cids);

      for (int i = Begin; i <= End; ++i) {
         const auto& s0i = shapeI0[i - Begin];
         const auto& s1i = shapeI1[i - Begin];
         const auto dsi = s1i - s0i;
         for (int j = Begin; j <= End; ++j) {
            const auto& s0j = shapeJ0[j - Begin];
            const auto& s1j = shapeJ1[j - Begin];
            const auto dsj = s1j - s0j;
            const auto tmp = -qA * (s0i * s0j + 0.5 * (dsi * s0j + s0i * dsj) + (1.0 / 3.0) * dsi * dsj);
            auto acc = 0.0;
            for (int k = Begin; k <= End - 1; ++k) {
               const auto& s0k = shapeK0[k - Begin];
               const auto& s1k = shapeK1[k - Begin];

               acc += (s1k - s0k) * tmp;
               const auto [x, y, z] = interp::rotateOrigin<D == 2 ? D : !D>(ci + i, cj + j, ck + k);
               #pragma omp atomic update
               J(x, y, z) += acc;
            } // end for(k)
         } // end for(j)
      } // end for(i)
   } // end deposit()


   template<int Support>
   static auto findRelayPoint(const auto& i0, const auto& i1, const auto& p1) {
      if constexpr (Support % 2 == 0) {
         return i0 == i1 ? p1 : std::max(i0, i1);
      } else {
         return i0 == i1 ? p1 : 0.5 * static_cast<double>(i1 + i0);
      }
   }


   static void updateJ(const auto& p, auto& emdata, const auto charge) {
      using Shape = interp::InterpolationShape<interpolation_order>::Type;
      static constexpr auto dtAxy = 1.0 / (dt * dx * dy);
      static constexpr auto dtAxz = 1.0 / (dt * dx * dz);
      static constexpr auto dtAyz = 1.0 / (dt * dy * dz);

      if (p.disabled) { return; }

      const auto x_coeff = p.weight * charge * dtAyz;
      const auto y_coeff = p.weight * charge * dtAxz;
      const auto z_coeff = p.weight * charge * dtAxy;

      const vec3<std::size_t> i0 = getCellIndices(p.old_location + 0.5);
      const vec3<std::size_t> i1 = getCellIndices(p.location + 0.5);

      const vec3<double> relay = {
         findRelayPoint<Shape::Support>(i0[0], i1[0], p.location[0]),
         findRelayPoint<Shape::Support>(i0[1], i1[1], p.location[1]),
         findRelayPoint<Shape::Support>(i0[2], i1[2], p.location[2]),
      };

      auto p0 = p.old_location - i0.as_type<double>();
      auto p1 = relay - i0.as_type<double>();

      auto s0i = Shape::shape_array(p0[0]);
      auto s0j = Shape::shape_array(p0[1]);
      auto s0k = Shape::shape_array(p0[2]);
      auto s1i = Shape::shape_array(p1[0]);
      auto s1j = Shape::shape_array(p1[1]);
      auto s1k = Shape::shape_array(p1[2]);

      deposit<0, Shape>(emdata.Jx, p0, p1, s0j, s0k, s0i, s1j, s1k, s1i, i0, x_coeff); // y-z-x
      deposit<1, Shape>(emdata.Jy, p0, p1, s0k, s0i, s0j, s1k, s1i, s1j, i0, y_coeff); // z-x-y
      deposit<2, Shape>(emdata.Jz, p0, p1, s0i, s0j, s0k, s1i, s1j, s1k, i0, z_coeff); // x-y-z

      if (i0 != i1) {
         p0 = relay - i1.as_type<double>();
         p1 = p.location - i1.as_type<double>();
         s0i = Shape::shape_array(p0[0]);
         s0j = Shape::shape_array(p0[1]);
         s0k = Shape::shape_array(p0[2]);
         s1i = Shape::shape_array(p1[0]);
         s1j = Shape::shape_array(p1[1]);
         s1k = Shape::shape_array(p1[2]);

         deposit<0, Shape>(emdata.Jx, p0, p1, s0j, s0k, s0i, s1j, s1k, s1i, i1, x_coeff); // y-z-x
         deposit<1, Shape>(emdata.Jy, p0, p1, s0k, s0i, s0j, s1k, s1i, s1j, i1, y_coeff); // z-x-y
         deposit<2, Shape>(emdata.Jz, p0, p1, s0i, s0j, s0k, s1i, s1j, s1k, i1, z_coeff); // x-y-z
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
