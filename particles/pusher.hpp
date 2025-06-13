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
   template<int D, int I>
   constexpr auto getShapeFunction(const auto x) {
      using CachedTSC = interp::Jit<interp::TSC>;
      using CachedCIC = interp::Jit<interp::CIC>;


      if constexpr ((D == 0 and I == 1) or (D == 1 and I != 2) or (D == 2 and I == 0)) {
         return CachedTSC{x};
      } else {
         return CachedCIC{x};
      }
   }
} // end namespace detail

template <int D, bool isB>
auto FieldToParticleInterp(const auto& F, auto p0, const auto& cids) {
   // static constexpr vec3 b0 = interp::rotateOrigin<D>(-1, (D == 1) != isB ? 0 : -1 , -1);
   // static constexpr vec3 b1 = interp::rotateOrigin<D>(D == 0 ? 0 : 1, 0, D == 2 ? 0 : 1);

   static constexpr vec3 b0{-1, -1 , -1};
   static constexpr vec3 b1 = interp::rotateOrigin<D>(
      D == 0 ? 0 : 1,
      (D == 1) == isB ? 0 : -1,
      D == 2 ? 0 : 1
   );

   // dbg(b0, b1);

   if constexpr (D == 0 or (D == 1 and isB)) {
      p0[0] -= 0.5;
   }
   if constexpr (D == 2 or (D == 1 and isB)) {
      p0[2] -= 0.5;
   }

   // // rotate elements to match the loop structure
   const auto& [ci, cj, ck] = interp::rotateOrigin<D>(cids);
   const auto& [x0, y0, z0] = interp::rotateOrigin<D>(p0);

   const auto shapeI = detail::getShapeFunction<D, 0>(x0);
   const auto shapeJ = detail::getShapeFunction<D, 1>(y0);
   const auto shapeK = detail::getShapeFunction<D, 2>(z0);

   // if constexpr (D == 2 and isB)
   //    dbg(shapeI.particle_position, shapeJ.particle_position, shapeK.particle_position);
   auto result = 0.0_fp;
   for (int i = b0[0]; i <= b1[0]; ++i) {
      const auto s0i = shapeI(i);
      for (int j = b0[1]; j <= b1[1]; ++j) {
         const auto s0j = shapeJ(j);
         for (int k = b0[2]; k <= b1[2]; ++k) {
            const auto s0k = shapeK(k);
            // undo the rotation to get proper indices back
            const auto [x, y, z] = interp::rotateOrigin<D == 2 ? D : !D>(ci + i, cj + j, ck + k);

            // dbg(x, y, z, F.dims());
            // if constexpr (D == 2 and isB)
            //    dbg(i, j, k, s0i, s0j, s0k, F(x, y, z));

            result += s0i * s0j * s0k * F(x, y, z);
            // if constexpr (D == 2 and isB)
            //    dbg(result);
            // dbg(result);
         } // end for(k)
      } // end for(j)
   } // end for(i)
   return result;
} // end FieldToParticle()


static std::array<double, 6> FieldAtParticle(Particle& p, const auto& emdata) {
   const auto cids = getCellIndices(p.location + 0.5);
   auto p0 = p.location - cids.as_type<compute_t>();
   // dbg(p.location, p0, cids);

   const auto exc = FieldToParticleInterp<0, 0>(emdata.Ex_total, p0, cids); // 1, 2, 0
   const auto eyc = FieldToParticleInterp<1, 0>(emdata.Ey_total, p0, cids); // 2, 0, 1
   const auto ezc = FieldToParticleInterp<2, 0>(emdata.Ez_total, p0, cids); // 0, 1, 2
   const auto bxc = FieldToParticleInterp<0, 1>(emdata.Bx_total, p0, cids); // 1, 2, 0
   const auto byc = FieldToParticleInterp<1, 1>(emdata.By_total, p0, cids); // 2, 0, 1
   const auto bzc = FieldToParticleInterp<2, 1>(emdata.Bz_total, p0, cids); // 0, 1, 2
   // dbg(exc, eyc, ezc, bxc, byc, bzc);
   return {exc, eyc, ezc, bxc, byc, bzc};
   // return {};
} // end FieldAtParticle

struct BorisPush {
   using emdata_t = electromagnetics::EMData;
   using group_t = ParticleGroup;

   static constexpr std::size_t BC_DEPTH = std::max(PMLDepth / 2zu, 1zu);

   static void update_velocity(Particle& p, const emdata_t& emdata, const auto qdt) {
      if (p.disabled) { return; }
      const auto emf = FieldAtParticle(p, emdata);

      const vec3<double> eps{qdt * emf[0], qdt * emf[1], qdt * emf[2]};
      const vec3<double> bet{qdt * emf[3], qdt * emf[4], qdt * emf[5]};

      const auto vm = p.velocity.as_type<double>() + eps;
      const auto t = bet / std::sqrt(1.0 + vm.length_squared() * constants::over_c_sqr<double>);
      const auto s = 2.0 * t / (1.0 + t.length_squared());
      const auto v = eps + (vm + cross(vm + cross(vm, t), s));
      p.gamma = std::sqrt(1.0 + v.length_squared() * constants::over_c_sqr<double>);
      p.velocity = v.as_type<compute_t>();
   } // end update_velocity()

   static void update_position(Particle& p) {
      static constexpr vec3 delta_inv{dt / dx, dt / dy, dt / dz};
      if (p.disabled) { return; }

      const auto new_loc = p.location + (delta_inv / p.gamma) * p.velocity;
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
