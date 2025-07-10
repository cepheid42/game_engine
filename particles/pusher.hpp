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
template <int D, typename Strategy>
auto FieldToParticleInterp(const auto& F, const auto& p0, const auto& cids) {
   using OuterShape = typename Strategy::OuterShape::type;
   using MiddleShape = typename Strategy::MiddleShape::type;
   using InnerShape = typename Strategy::InnerShape::type;

   static constexpr auto iBegin = OuterShape::Begin;
   static constexpr auto jBegin = MiddleShape::Begin;
   static constexpr auto kBegin = InnerShape::Begin;
   static constexpr auto iEnd = OuterShape::End;
   static constexpr auto jEnd = MiddleShape::End;
   static constexpr auto kEnd = InnerShape::End;

   const auto& [ci, cj, ck] = interp::rotateOrigin<D>(cids);
   const auto& [x0, y0, z0] = interp::rotateOrigin<D>(p0);

   const auto shapeI = OuterShape::shape_array(x0);
   const auto shapeJ = MiddleShape::shape_array(y0);
   const auto shapeK = InnerShape::shape_array(z0);

   auto result = 0.0_fp;
   for (int i = iBegin; i <= iEnd; ++i) {
      const auto& s0i = shapeI[i - iBegin];
      for (int j = jBegin; j <= jEnd; ++j) {
         const auto& s0j = shapeJ[j - jBegin];
         for (int k = kBegin; k <= kEnd - 1; ++k) {
            const auto& s0k = shapeK[k - kBegin];
            const auto [x, y, z] = interp::rotateOrigin<D == 2 ? D : !D>(ci + i, cj + j, ck + k);
            result += s0i * s0j * s0k * F(x, y, z);
         } // end for(k)
      } // end for(j)
   } // end for(i)
   return result;
} // end FieldToParticle()

static std::array<double, 6> fieldAtParticle(const Particle& p, const auto& emdata) {
   static_assert(interpolation_order != 0);
   using AssignmentShape = interp::InterpolationShape<interpolation_order>;
   using ReducedShape = interp::InterpolationShape<interpolation_order - 1>;
   using EStrategy = interp::InterpolationStrategy<AssignmentShape, AssignmentShape, ReducedShape>; // 1, 2, 0
   using BStrategy = interp::InterpolationStrategy<ReducedShape, ReducedShape, AssignmentShape>; // 1, 2, 0

   const auto cids = getCellIndices(p.location);
   auto p0 = p.location - cids.as_type<compute_t>();

   const auto exc = FieldToParticleInterp<0, EStrategy>(emdata.Ex_total, p0, cids);
   const auto eyc = FieldToParticleInterp<1, EStrategy>(emdata.Ey_total, p0, cids);
   const auto ezc = FieldToParticleInterp<2, EStrategy>(emdata.Ez_total, p0, cids);
   const auto bxc = FieldToParticleInterp<0, BStrategy>(emdata.Bx_total, p0, cids);
   const auto byc = FieldToParticleInterp<1, BStrategy>(emdata.By_total, p0, cids);
   const auto bzc = FieldToParticleInterp<2, BStrategy>(emdata.Bz_total, p0, cids);

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

      // u = gamma * v
      const auto um = p.gamma * p.velocity + eps;
      const auto t = bet / std::sqrt(1.0 + um.length_squared() * constants::over_c_sqr<double>);
      const auto s = 2.0 * t / (1.0 + t.length_squared());
      const auto u_prime = um + cross(um, t);
      const auto u_plus = um + cross(u_prime, s);
      const auto u = u_plus + eps;

      p.gamma = calculateGammaV(u);
      p.velocity = u / p.gamma;
   } // end update_velocity()

   static void update_position(Particle& p) {
      // todo: If I store v/c, these could be dt*c/dx = cfl, which would be easy...
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
