#ifndef PUSHER_HPP
#define PUSHER_HPP

#include "program_params.hpp"
#include "particles.hpp"
#include "em_data.hpp"
#include "constants.hpp"
#include "interpolation.hpp"

// #include "dbg.h"

#include <cmath>
#include <cassert>

namespace tf::particles {
namespace detail {
   template<int C, bool isB>
   constexpr auto getShapeFunction(const auto x) {
      using NGP = interp::Jit<interp::NGP>;
      using CIC = interp::Jit<interp::CIC>;
      // using TSC = interp::Jit<interp::TSC>;

      if constexpr (C == 2) {
         if constexpr (isB) {
            return CIC{x};
         } else {
            return NGP{x};
         }
      } else {
         if constexpr (isB) {
            return NGP{x};
         } else {
            return CIC{x};
         }
      }
   }
} // end namespace detail

template <int D, bool isB>
auto FieldToParticleInterp(const auto& F, auto p0, const auto& cids, const bool print=false) {
   static constexpr vec3 b0{0, 0, 0};
   // Ex/Ey/Ez -> {1, 1, 0}, Bx/By/Bz -> {0, 0, 1}
   static constexpr vec3 b1{
      isB ? 0 : 1,
      isB ? 0 : 1,
      isB ? 1 : 0
   };

   // static constexpr vec3 b0{-1, -1, -1};
   // // Ex/Ey/Ez -> {1, 1, 0}, Bx/By/Bz -> {0, 0, 1}
   // static constexpr vec3 b1{
   //    isB ? -1 : 0,
   //    isB ? -1 : 0,
   //    isB ? 0 : -1
   // };

   // // rotate elements to match the loop structure
   const auto& [ci, cj, ck] = interp::rotateOrigin<D>(cids);
   const auto& [x0, y0, z0] = interp::rotateOrigin<D>(p0);

   const auto shapeI = detail::getShapeFunction<0, isB>(x0);
   const auto shapeJ = detail::getShapeFunction<1, isB>(y0);
   const auto shapeK = detail::getShapeFunction<2, isB>(z0);

   auto result = 0.0_fp;
   for (int i = b0[0]; i <= b1[0]; ++i) {
      const auto s0i = shapeI(i);
      if (print)
         std::println("i = {}: {}", i, s0i);

      for (int j = b0[1]; j <= b1[1]; ++j) {
         const auto s0j = shapeJ(j);
         if (print)
            std::println("j = {}: {}", j, s0j);

         for (int k = b0[2]; k <= b1[2]; ++k) {
            const auto s0k = shapeK(k);
            if (print)
               std::println("k = {}: {}", k, s0k);

            // undo the rotation to get proper indices back
            const auto [x, y, z] = interp::rotateOrigin<D == 2 ? D : !D>(ci + i, cj + j, ck + k);
            result += s0i * s0j * s0k * F(x, y, z);
            if (print)
               std::println("{}, {}, {}: {}", x, y, z, result);
         } // end for(k)
      } // end for(j)
   } // end for(i)

   return result;
} // end FieldToParticle()


static std::array<double, 6> fieldAtParticle(const Particle& p, const auto& emdata) {
   const auto cids = getCellIndices(p.location);
   // const auto cids = getCellIndices(p.location + 1.0);
   auto p0 = p.location - cids.as_type<compute_t>();

   // std::println("{}, {}", p0, cids);

   // std::println("Ex");
   const auto exc = FieldToParticleInterp<0, 0>(emdata.Ex_total, p0, cids); // 1, 2, 0
   // std::println("Ey");
   const auto eyc = FieldToParticleInterp<1, 0>(emdata.Ey_total, p0, cids); // 2, 0, 1
   // std::println("Ez");
   const auto ezc = FieldToParticleInterp<2, 0>(emdata.Ez_total, p0, cids); // 0, 1, 2
   // std::println("Bx");
   const auto bxc = FieldToParticleInterp<0, 1>(emdata.Bx_total, p0, cids); // 1, 2, 0
   // std::println("By");
   const auto byc = FieldToParticleInterp<1, 1>(emdata.By_total, p0, cids); // 2, 0, 1
   // std::println("Bz");
   const auto bzc = FieldToParticleInterp<2, 1>(emdata.Bz_total, p0, cids); // 0, 1, 2
   // return {};
   return {exc, eyc, ezc, bxc, byc, bzc};
} // end FieldAtParticle

struct BorisPush {
   using emdata_t = electromagnetics::EMData;
   using group_t = ParticleGroup;

   static constexpr std::size_t BC_DEPTH = std::max(PMLDepth / 2zu, 1zu);

   static void update_velocity(Particle& p, const emdata_t& emdata, const auto qdt) {
      if (p.disabled) { return; }
      const auto emf = fieldAtParticle(p, emdata);

      const vec3<double> eps{qdt * emf[0], qdt * emf[1], qdt * emf[2]};
      const vec3<double> bet{qdt * emf[3], qdt * emf[4], qdt * emf[5]};

      const auto vm = p.gamma * p.velocity.as_type<double>() + eps;
      const auto t = bet / std::sqrt(1.0 + vm.length_squared() * constants::over_c_sqr<double>);
      const auto s = 2.0 * t / (1.0 + t.length_squared());
      const auto v = eps + (vm + cross(vm + cross(vm, t), s));
      p.gamma = std::sqrt(1.0 + v.length_squared() * constants::over_c_sqr<double>);
      p.velocity = v.as_type<compute_t>() / p.gamma;
   } // end update_velocity()

   static void update_position(Particle& p) {
      static constexpr vec3 delta_inv{dt / dx, dt / dy, dt / dz};
      if (p.disabled) { return; }

      const auto new_loc = p.location + delta_inv * p.velocity;
      const auto [inew, jnew, knew] = getCellIndices(new_loc);

      p.disabled = inew < BC_DEPTH or inew > Nx - BC_DEPTH or knew < BC_DEPTH or knew > Nz - BC_DEPTH;
      p.old_location = p.location;
      p.location = new_loc;
   } // end update_position()


   static void advance_velocity(group_t& g, const emdata_t& emdata) {
      #pragma omp parallel for num_threads(nThreads)
      for (std::size_t pid = 0; pid < g.num_particles(); pid++) {
         update_velocity(g.particles[pid], emdata, g.qdt_over_2m);
      }
   } // end advance_velocity

   static void advance_position(group_t& g) {
      #pragma omp parallel for num_threads(nThreads)
      for (std::size_t pid = 0; pid < g.num_particles(); pid++) {
         update_position(g.particles[pid]);
      }
   } // end advance_position

   static void backstep_velocity(group_t& g, const emdata_t& emdata) {
      #pragma omp parallel for num_threads(nThreads)
      for (std::size_t pid = 0; pid < g.num_particles(); pid++) {
         update_velocity(g.particles[pid], emdata, -0.5 * g.qdt_over_2m);
      }
   }

   // static void backstep_position(group_t& g) {
   //    static constexpr vec3 delta_inv{dt / dx, dt / dy, dt / dz};
   //
   //    #pragma omp parallel for num_threads(nThreads)
   //    for (std::size_t pid = 0; pid < g.num_particles(); pid++) {
   //       auto& p = g.particles[pid];
   //       const auto old_loc = p.location + (-0.5 * delta_inv / p.gamma) * p.velocity;
   //       const auto [inew, jnew, knew] = getCellIndices(old_loc);
   //       p.disabled = inew < BC_DEPTH or inew > Nx - BC_DEPTH or knew < BC_DEPTH or knew > Nz - BC_DEPTH;
   //       p.old_location = old_loc;
   //    }
   // }

   static void advance(group_t& g, const emdata_t& emdata, const std::size_t step) {
      advance_velocity(g, emdata);
      advance_position(g);
      if (step % group_t::SORT_INTERVAL == 0) {
         g.sort_particles();
      }
   } // end operator()
}; // end struct BorisPush
} // end namespace tf::particles


#endif //PUSHER_HPP
