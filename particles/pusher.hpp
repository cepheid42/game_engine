#ifndef PUSHER_HPP
#define PUSHER_HPP

#include "program_params.hpp"
#include "particles.hpp"
#include "em_data.hpp"
#include "constants.hpp"
#include "interpolation.hpp"

// #include "dbg.h"

#include <cmath>

namespace tf::particles {
template <int D, bool isB>
auto FieldToParticleInterp(const auto& F, const auto& cids,
                           const auto& shapeI, const auto& shapeJ, const auto& shapeK) {
   static constexpr vec3 b0{-1, -1, -1};
   static constexpr vec3 b1 = interp::rotateOrigin<D>(1, (D == 1) != isB ? -1 : 0, 1);

   // // rotate elements to match the loop structure
   const auto [ci, cj, ck] = interp::rotateOrigin<D>(cids);

   auto result = 0.0_fp;
   for (int i = b0[0]; i <= b1[0]; ++i) {
      const auto s0i = shapeI(i);
      for (int j = b0[1]; j <= b1[1]; ++j) {
         const auto s0j = shapeJ(j);
         for (int k = b0[2]; k <= b1[2]; ++k) {
            const auto s0k = shapeK(k);
            // undo the rotation to get proper indices back
            const auto [x, y, z] = interp::rotateOrigin<D == 2 ? D : !D>(i + ci, j + cj, k + ck);
            result += s0i * s0j * s0k * F(x, y, z);
         } // end for(k)
      } // end for(j)
   } // end for(i)
   return result;
} // end FieldToParticle()


static std::array<double, 6> FieldAtParticle(Particle& p, const auto& emdata) {
   using TSCCache = interp::Jit<interp::TSC>;
   using CICCache = interp::Jit<interp::CIC>;

   // const auto new_loc = p.location + 0.5;
   // const auto cids = getCellIndices(new_loc);
   // const auto p0 = new_loc - cids.as_type<compute_t>();
   // // todo: does p0[1] need to be corrected?

   const auto cids = getCellIndices(p.location);
   auto p0 = p.location - cids.as_type<compute_t>();
   p0[1] -= 0.5_fp;

   const TSCCache shapeI(p0[0]);
   const CICCache shapeJ(p0[1]);
   const TSCCache shapeK(p0[2]);

   const auto exc = FieldToParticleInterp<0, 0>(emdata.Ex_total, cids, shapeJ, shapeK, shapeI); // 1, 2, 0
   const auto eyc = FieldToParticleInterp<1, 0>(emdata.Ey_total, cids, shapeK, shapeI, shapeJ); // 2, 0, 1
   const auto ezc = FieldToParticleInterp<2, 0>(emdata.Ez_total, cids, shapeI, shapeJ, shapeK); // 0, 1, 2
   const auto bxc = FieldToParticleInterp<0, 1>(emdata.Bx_total, cids, shapeJ, shapeK, shapeI); // 1, 2, 0
   const auto byc = FieldToParticleInterp<1, 1>(emdata.By_total, cids, shapeK, shapeI, shapeJ); // 2, 0, 1
   const auto bzc = FieldToParticleInterp<2, 1>(emdata.Bz_total, cids, shapeI, shapeJ, shapeK); // 0, 1, 2

   return {exc, eyc, ezc, bxc, byc, bzc};
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

      // const auto vminus = p.velocity + eps;
      // const auto t = bet / std::sqrt(1.0 + vminus.length_squared() * constants::over_c_sqr<double>);
      // const auto s = 2.0 * t / (1.0 + t.length_squared());
      // const auto v_prime = vminus + cross(vminus, t);
      // const auto vplus = vminus + cross(v_prime, s);
      // p.velocity = vplus + eps;
      // p.gamma = std::sqrt(1.0 + p.velocity.length_squared() * constants::over_c_sqr<double>);
   } // end update_velocity()

   static void update_position(Particle& p) {
      static constexpr vec3 delta_inv{dt / dx, dt / dy, dt / dz};
      if (p.disabled) { return; }

      const auto new_loc = p.location + delta_inv * p.velocity / p.gamma;
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

   static void operator()(group_t& g, const emdata_t& emdata, const std::size_t step) {
      advance_velocity(g, emdata);
      advance_position(g);

      if (step % group_t::SORT_INTERVAL == 0) {
         g.sort_particles();
      }
   } // end operator()
}; // end struct BorisPush
} // end namespace tf::particles


#endif //PUSHER_HPP
