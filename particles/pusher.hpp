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
   using OShape = typename Strategy::OuterShape::type;
   using MShape = typename Strategy::MiddleShape::type;
   using IShape = typename Strategy::InnerShape::type;

   static constexpr auto iBegin = OShape::Begin;
   static constexpr auto jBegin = MShape::Begin;
   static constexpr auto kBegin = IShape::Begin;
   static constexpr auto iEnd = OShape::End;
   static constexpr auto jEnd = MShape::End;
   static constexpr auto kEnd = IShape::End;

   const auto& [ci, cj, ck] = interp::rotateOrigin<D>(cids);
   const auto& [x0, y0, z0] = interp::rotateOrigin<D>(p0);

   const auto shapeI = OShape::shape_array(x0);
   const auto shapeJ = MShape::shape_array(y0);
   const auto shapeK = IShape::shape_array(z0);

   auto avg = 0.0;

   auto result = 0.0_fp;
   for (int i = iBegin; i <= iEnd; ++i) {
      const auto& s0i = shapeI[i - iBegin];
      for (int j = jBegin; j <= jEnd; ++j) {
         const auto& s0j = shapeJ[j - jBegin];
         for (int k = kBegin; k <= kEnd; ++k) {
            const auto& s0k = shapeK[k - kBegin];
            avg += s0i * s0j * s0k;
            const auto& [x, y, z] = interp::rotateOrigin<D == 2 ? D : !D>(ci + i, cj + j, ck + k);
            result += s0i * s0j * s0k * F(x, y, z);
         } // end for(k)
      } // end for(j)
   } // end for(i)

   return result;
} // end FieldToParticle()

static std::array<vec3<double>, 2> fieldAtParticle(const Particle& p, const auto& emdata, const auto qdt) {
   using AssignmentShape = interp::InterpolationShape<interpolation_order>;
   using ReducedShape = interp::InterpolationShape<interpolation_order - 1>;
   using EStrategy = interp::InterpolationStrategy<AssignmentShape, AssignmentShape, ReducedShape>;
   using BStrategy = interp::InterpolationStrategy<ReducedShape, ReducedShape, AssignmentShape>;

   const vec3 cid = getCellIndices<compute_t>(p.location + 0.5);

   // Ex H,F,F
   const vec3 ExNloc{cid[0] - 0.5, cid[1], cid[2]};
   const auto ExP0 = p.location - ExNloc;
   const auto exc = FieldToParticleInterp<0, EStrategy>(emdata.Ex_total, ExP0, ExNloc.as_type<std::size_t>());

   // Ey F,H,F
   const vec3 EyNloc{cid[0], cid[1] - 0.5, cid[2]};
   const auto EyP0 = p.location - EyNloc;
   const auto eyc = FieldToParticleInterp<1, EStrategy>(emdata.Ey_total, EyP0, EyNloc.as_type<std::size_t>());

   // Ez
   const vec3 EzNloc{cid[0], cid[1], cid[2] - 0.5};
   const auto EzP0 = p.location - EzNloc;
   const auto ezc = FieldToParticleInterp<2, EStrategy>(emdata.Ez_total, EzP0, EzNloc.as_type<std::size_t>());

   // Bx
   const vec3 BxNloc{cid[0], cid[1] - 0.5, cid[2] - 0.5};
   const auto BxP0 = p.location - BxNloc;
   const auto bxc = FieldToParticleInterp<0, BStrategy>(emdata.Bx_total, BxP0, BxNloc.as_type<std::size_t>());

   // By
   const vec3 ByNloc{cid[0] - 0.5, cid[1], cid[2] - 0.5};
   const auto ByP0 = p.location - ByNloc;
   const auto byc = FieldToParticleInterp<1, BStrategy>(emdata.By_total, ByP0, ByNloc.as_type<std::size_t>());

   // Bz
   const vec3 BzNloc{cid[0] - 0.5, cid[1] - 0.5, cid[2]};
   const auto BzP0 = p.location - BzNloc;
   const auto bzc = FieldToParticleInterp<2, BStrategy>(emdata.Bz_total, BzP0, BzNloc.as_type<std::size_t>());

   return {qdt * vec3{exc, eyc, ezc}, qdt * vec3{bxc, byc, bzc}};
} // end FieldAtParticle

struct BorisPush {
   using emdata_t = electromagnetics::EMData;
   using group_t = ParticleGroup;

   static constexpr std::size_t BC_DEPTH = std::max(PMLDepth / 2zu, 1zu);

   static void update_velocity(Particle& p, const emdata_t& emdata, const auto qdt) {
      if (p.disabled) { return; }
      const auto& [eps, bet] = fieldAtParticle(p, emdata, qdt);

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
