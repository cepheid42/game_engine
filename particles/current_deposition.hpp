#ifndef CURRENT_DEPOSITION_HPP
#define CURRENT_DEPOSITION_HPP

#include "interpolation.hpp"
#include "particles.hpp"

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


   template<int Order>
   static auto findRelayPoint(const auto i0, const auto i1, const auto p1) {
      if constexpr (Order % 2 == 0) {
         return i0 == i1 ? p1 : 0.5 * static_cast<double>(i1 + i0);
      } else {
         return i0 == i1 ? p1 : std::fmax(i0, i1);
      }
   }


   static void updateJ(const auto& p, auto& emdata, const auto charge) {
      using Shape = interp::InterpolationShape<interpolation_order>::Type;

      if (p.disabled) { return; }

      const auto x_coeff = p.weight * charge / (dt * dy * dz);
      const auto y_coeff = p.weight * charge / (dt * dx * dz);
      const auto z_coeff = p.weight * charge / (dt * dx * dz);

      constexpr auto offset = Shape::Order % 2 == 0 ? 0.5 : 0.0;
      const vec3<std::size_t> i0 = getCellIndices(p.old_location + offset);
      const vec3<std::size_t> i1 = getCellIndices(p.location + offset);

      const vec3<double> relay = {
         findRelayPoint<Shape::Order>(i0[0], i1[0], p.location[0]),
         findRelayPoint<Shape::Order>(i0[1], i1[1], p.location[1]),
         findRelayPoint<Shape::Order>(i0[2], i1[2], p.location[2]),
      };

      auto p0 = p.old_location - i0.as_type<double>();
      auto p1 = relay - i0.as_type<double>();

      auto s0i = Shape::shape_array(p0[0]);
      auto s0j = Shape::shape_array(p0[1]);
      auto s0k = Shape::shape_array(p0[2]);
      auto dsi = Shape::ds_array(p1[0], s0i);
      auto dsj = Shape::ds_array(p1[1], s0j);
      auto dsk = Shape::ds_array(p1[2], s0k);

      // std::println("{}, {}", p.old_location, p.location);
      // std::println("{}, {}, {}", i0, i1, relay);
      // std::println("{}, {}", p0, p1);
      // auto s1i = Shape::shape_array(p1[0]);
      // auto s1j = Shape::shape_array(p1[1]);
      // auto s1k = Shape::shape_array(p1[2]);
      // // First Order
      // std::println("s0i: ({}, {})", s0i[0], s0i[1]);
      // std::println("s0j: ({}, {})", s0j[0], s0j[1]);
      // std::println("s0k: ({}, {})", s0k[0], s0k[1]);
      // std::println("s1i: ({}, {})", s1i[0], s1i[1]);
      // std::println("s1j: ({}, {})", s1j[0], s1j[1]);
      // std::println("s1k: ({}, {})", s1k[0], s1k[1]);

      // // Second Order
      // std::println("s0i: ({}, {}, {})", s0i[0], s0i[1], s0i[2]);
      // std::println("s0j: ({}, {}, {})", s0j[0], s0j[1], s0j[2]);
      // std::println("s0k: ({}, {}, {})", s0k[0], s0k[1], s0k[2]);
      // std::println("s1i: ({}, {}, {})", s1i[0], s1i[1], s1i[2]);
      // std::println("s1j: ({}, {}, {})", s1j[0], s1j[1], s1j[2]);
      // std::println("s1k: ({}, {}, {})", s1k[0], s1k[1], s1k[2]);

      deposit<0, Shape>(emdata.Jx, p0[0], p1[0], s0j, s0k, dsj, dsk, dsi, i0[1], i0[2], i0[0], x_coeff); // y-z-x
      deposit<1, Shape>(emdata.Jy, p0[1], p1[1], s0k, s0i, dsk, dsi, dsj, i0[2], i0[0], i0[1], y_coeff); // z-x-y
      deposit<2, Shape>(emdata.Jz, p0[2], p1[2], s0i, s0j, dsi, dsj, dsk, i0[0], i0[1], i0[2], z_coeff); // x-y-z

      if (i0 != i1) {
         p0 = relay - i1.as_type<double>();
         p1 = p.location - i1.as_type<double>();
         s0i = Shape::shape_array(p0[0]);
         s0j = Shape::shape_array(p0[1]);
         s0k = Shape::shape_array(p0[2]);
         dsi = Shape::ds_array(p1[0], s0i);
         dsj = Shape::ds_array(p1[1], s0j);
         dsk = Shape::ds_array(p1[2], s0k);

         // std::println("Cell crossing!");
         // std::println("{}, {}", p0, p1);
         // s1i = Shape::shape_array(p1[0]);
         // s1j = Shape::shape_array(p1[1]);
         // s1k = Shape::shape_array(p1[2]);
         // // First Order
         // std::println("s0i: ({}, {})", s0i[0], s0i[1]);
         // std::println("s0j: ({}, {})", s0j[0], s0j[1]);
         // std::println("s0k: ({}, {})", s0k[0], s0k[1]);
         // std::println("s1i: ({}, {})", s1i[0], s1i[1]);
         // std::println("s1j: ({}, {})", s1j[0], s1j[1]);
         // std::println("s1k: ({}, {})", s1k[0], s1k[1]);
         // // Second Order
         // std::println("s0i: ({}, {}, {})", s0i[0], s0i[1], s0i[2]);
         // std::println("s0j: ({}, {}, {})", s0j[0], s0j[1], s0j[2]);
         // std::println("s0k: ({}, {}, {})", s0k[0], s0k[1], s0k[2]);
         // std::println("s1i: ({}, {}, {})", s1i[0], s1i[1], s1i[2]);
         // std::println("s1j: ({}, {}, {})", s1j[0], s1j[1], s1j[2]);
         // std::println("s1k: ({}, {}, {})", s1k[0], s1k[1], s1k[2]);

         deposit<0, Shape>(emdata.Jx, p0[0], p1[0], s0j, s0k, dsj, dsk, dsi, i1[1], i1[2], i1[0], x_coeff); // y-z-x
         deposit<1, Shape>(emdata.Jy, p0[1], p1[1], s0k, s0i, dsk, dsi, dsj, i1[2], i1[0], i1[1], y_coeff); // z-x-y
         deposit<2, Shape>(emdata.Jz, p0[2], p1[2], s0i, s0j, dsi, dsj, dsk, i1[0], i1[1], i1[2], z_coeff); // x-y-z
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
