#ifndef PUSHER_HPP
#define PUSHER_HPP

#include "program_params.hpp"
#include "em_params.hpp"
#include "particles.hpp"
#include "em_data.hpp"
#include "constants.hpp"
#include "interpolation.hpp"

#include "dbg.h"

#include <cmath>

namespace tf::particles
{
template<int D, bool isB>
auto FieldToParticleInterp(const auto& F, const auto& p0, const auto& cids)
{
   using CachedShape = interp::Jit<interp::TSC>;
   // (x, y, z) -> (1, 2, 0) -> (y, z, x) | Ex/Bx D = 0
   // (x, y, z) -> (2, 0, 1) -> (z, x, y) | Ey/By D = 1
   // (x, y, z) ->           -> (x, y, z) | Ez D = 2

   static constexpr vec3 b0 = interp::rotateOrigin<D>(vec3{-1, 0, -1});
   static constexpr vec3 b1 = interp::rotateOrigin<D>(vec3{1, static_cast<int>((D == 1) == isB), 1});
   // b1[1] is 0 for Ey/Bx/Bz and 1 for Ex/Ez/By

   // rotate elements to match loop structure
   const auto [ci, cj, ck] = interp::rotateOrigin<D>(cids);
   const auto [x0, y0, z0] = interp::rotateOrigin<D>(p0);

   const CachedShape shapeI(x0);
   const CachedShape shapeJ(y0);
   const CachedShape shapeK(z0);

   auto result = 0.0_fp;
   for (int i = b0[0]; i <= b1[0]; ++i)
   {
      const auto s0i = shapeI(ci + i);

      for (int j = b0[1]; j <= b1[1]; ++j)
      {
         const auto s0j = shapeJ(cj + j);

         for (int k = b0[2]; k <= b1[2]; ++k)
         {
            const auto s0k = shapeK(ck + k);

            // undo the rotation to get proper indices back
            const auto [x, y, z] = interp::rotateOrigin<D == 2 ? D : !D>(vec3{i + ci, j + cj, k + ck});
            result += s0i * s0j * s0k * F(x, y, z);
            // dbg(s0i, s0j, s0k, F(x, y, z));
         } // end for(k)
      }    // end for(j)
   }       // end for(i)
   return result;
} // end EFieldToParticle()

static std::array<double, 6> FieldAtParticle(Particle& p, const auto& emdata)
{
   const auto cids = getCIDs(p.location);

   assertm(cids[0] < Ncx - 1 and cids[0] > 1, "FieldAtParticle: Invalid particle cell in x.");
   // assertm(cids[1] < Ncy - 1 and cids[1] > 1, "FieldAtParticle: Invalid particle cell in y.");
   assertm(cids[2] < Ncz - 1 and cids[2] > 1, "FieldAtParticle: Invalid particle cell in z.");

   /*
    *
    * todo: There has to be a way to create the interpolation stencils on a per-cell basis
    *       so I can save time making the shapes for each and every particle.
    *
    */
   const auto exc = FieldToParticleInterp<0, false>(emdata.Ex, p.location, cids);
   const auto eyc = FieldToParticleInterp<1, false>(emdata.Ey, p.location, cids);
   const auto ezc = FieldToParticleInterp<2, false>(emdata.Ez, p.location, cids);

   const auto bxc = FieldToParticleInterp<0, true>(emdata.Bx, p.location, cids);
   const auto byc = FieldToParticleInterp<1, true>(emdata.By, p.location, cids);
   const auto bzc = FieldToParticleInterp<2, true>(emdata.Bz, p.location, cids);

   return {exc, eyc, ezc, bxc, byc, bzc};
} // end FieldAtParticle

struct BorisPush
{
   using emdata_t = electromagnetics::EMData;
   using group_t  = ParticleGroup;

   static constexpr std::size_t BC_DEPTH = PMLDepth / 2;

   static void update_velocity(Particle& p, const emdata_t& emdata, auto qdt)
   {
      if (p.disabled) { return; }
      const auto emf = FieldAtParticle(p, emdata);

      // dbg(emf);

      const vec3<double> eps{qdt * emf[0], qdt * emf[1], qdt * emf[2]};
      const vec3<double> bet{qdt * emf[3], qdt * emf[4], qdt * emf[5]};

      const auto vm = p.velocity.as_type<double>() + eps;
      const auto t  = bet / std::sqrt(1.0 + vm.length_squared() * constants::over_c_sqr<double>);
      const auto s  = 2.0 * t / (1.0 + t.length_squared());
      const auto v  = eps + (vm + cross(vm + cross(vm, t), s));
      p.gamma       = std::sqrt(1.0 + v.length_squared() * constants::over_c_sqr<double>);
      p.velocity    = v.as_type<compute_t>();
   } // end update_velocity()

   // #pragma omp declare simd notinbranch
   static void update_position(Particle& p)
   {
      static constexpr vec3 delta_inv{dt / dx, dt / dy, dt / dz};
      if (p.disabled) { return; }

      const auto new_loc            = p.location + delta_inv * p.velocity / p.gamma;
      const auto [inew, jnew, knew] = getCIDs(new_loc);

      p.disabled     = inew < BC_DEPTH or inew > Nx - BC_DEPTH or knew < BC_DEPTH or knew > Nz - BC_DEPTH;
      p.old_location = p.location;
      p.location     = new_loc;
      // dbg(p.old_location, p.location);
   } // end update_position()


   static void advance_velocity(group_t& g, const emdata_t& emdata)
   {
      #pragma omp parallel for num_threads(nThreads)
      for (std::size_t pid = 0; pid < g.num_particles(); pid++)
      {
         update_velocity(g.particles[pid], emdata, g.qdt_over_2m);
      }
   } // end advance_velocity

   static void advance_position(group_t& g)
   {
      #pragma omp parallel for simd num_threads(nThreads)
      for (std::size_t pid = 0; pid < g.num_particles(); pid++)
      {
         update_position(g.particles[pid]);
      }
      // std::erase_if(g.particles, [](const Particle& p) { return p.disabled; });
   } // end advance_position


   static void backstep_velocity(group_t& g, const emdata_t& emdata)
   {
      // Easiest way is to just copy the velocity update and make qdt negative
      #pragma omp parallel for num_threads(nThreads)
      for (std::size_t pid = 0; pid < g.num_particles(); pid++)
      {
         update_velocity(g.particles[pid], emdata, -0.5 * g.qdt_over_2m);
      }
   }

   static void operator()(group_t& g, const emdata_t& emdata, const std::size_t step)
   {
      advance_velocity(g, emdata);
      advance_position(g);

      if (step % group_t::SORT_INTERVAL == 0)
      {
         g.sort_particles();
      }
   } // end operator()
};   // end struct BorisPush
}    // end namespace tf::particles


#endif //PUSHER_HPP
