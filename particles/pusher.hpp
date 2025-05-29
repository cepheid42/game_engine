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

namespace tf::particles {
   auto ExInterp(const auto& Ex, const auto& p0, const auto& cids) {
      using jitTSC = interp::Jit<interp::TSC>;

      const auto [x0, y0, z0] = p0;
      const auto [ci, cj, ck] = cids;
      const jitTSC shapeI(x0 - ci);
      const jitTSC shapeJ(y0 - cj);
      const jitTSC shapeK(z0 - ck);

      auto result = 0.0_fp;
      for (int i = -1; i <= 1; ++i) {
         const auto s0i = shapeI(i);
         for (int j = -1; j <= 0; ++j) {
            const auto s0j = shapeJ(j);
            for (int k = -1; k <= 1; ++k) {
               const auto s0k = shapeK(k);
               result += s0i * s0j * s0k * Ex(ci + i, cj + j, ck + k);
            } // end for(k)
         } // end for(j)
      } // end for(i)
      return result;
   }

   auto EyInterp(const auto& Ey, const auto& p0, const auto& cids) {
      using jitTSC = interp::Jit<interp::TSC>;
      using jitCIC = interp::Jit<interp::CIC>;

      const auto [x0, y0, z0] = p0;
      const auto [ci, cj, ck] = cids;
      const jitTSC shapeI(x0 - ci);
      const jitCIC shapeJ(y0 - cj);
      const jitTSC shapeK(z0 - ck);

      auto result = 0.0_fp;
      for (int i = -1; i <= 1; ++i) {
         const auto s0i = shapeI(i);
         for (int j = -1; j <= -1; ++j) {
            const auto s0j = shapeJ(j);
            for (int k = -1; k <= 1; ++k) {
               const auto s0k = shapeK(k);
               result += s0i * s0j * s0k * Ey(ci + i, cj + j, ck + k);
            } // end for(k)
         } // end for(j)
      } // end for(i)

      return result;
   }

   auto EzInterp(const auto& Ez, const auto& p0, const auto& cids) {
      using cachedTSC = interp::Cached<interp::TSC>;
      using cachedCIC = interp::Cached<interp::CIC>;

      const auto [x0, y0, z0] = p0;
      const auto [ci, cj, ck] = cids;

      const cachedTSC shapeI(x0 - ci);
      const cachedCIC shapeJ(y0 - cj);
      const cachedTSC shapeK(z0 - ck);

      auto result = 0.0_fp;
      for (int i = -1; i <= 1; ++i) {
         const auto s0i = shapeI(i);
         for (int j = -1; j <= 0; ++j) {
            const auto s0j = shapeJ(j);
            for (int k = -1; k <= 1; ++k) {
               const auto s0k = shapeK(k);
               result += s0i * s0j * s0k * Ez(ci + i, cj + j, ck + k);
            } // end for(k)
         } // end for(j)
      } // end for(i)
      return result;
   }

   auto BxInterp(const auto& Bx, const auto& p0, const auto& cids) {
      using jitTSC = interp::Jit<interp::TSC>;
      using jitCIC = interp::Jit<interp::CIC>;

      const auto [x0, y0, z0] = p0;
      const auto [ci, cj, ck] = cids;
      const jitTSC shapeI(x0 - ci);
      const jitCIC shapeJ(y0 - cj);
      const jitTSC shapeK(z0 - ck);

      auto result = 0.0_fp;
      for (int i = -1; i <= 1; ++i) {
         const auto s0i = shapeI(i);
         for (int j = -1; j <= -1; ++j) {
            const auto s0j = shapeJ(j);
            for (int k = -1; k <= 1; ++k) {
               const auto s0k = shapeK(k);
               result += s0i * s0j * s0k * Bx(ci + i, cj + j, ck + k);
            } // end for(k)
         } // end for(j)
      } // end for(i)

      // assert(result == tmp);
      return result;
   }

   auto ByInterp(const auto& By, const auto& p0, const auto& cids) {
      using jitTSC = interp::Jit<interp::TSC>;
      using jitCIC = interp::Jit<interp::CIC>;
      using cachedTSC = interp::Cached<interp::TSC>;

      const auto [x0, y0, z0] = p0;
      const auto [ci, cj, ck] = cids;

      const jitTSC shapeI(x0 - ci);
      const jitCIC shapeJ(y0 - cj);
      const jitTSC shapeK(z0 - ck);

      auto result = 0.0_fp;
      for (int i = -1; i <= 1; ++i) {
         const auto s0i = shapeI(i);
         // assert(s0i == rshapeI(i));

         for (int j = -1; j <= 0; ++j) {
            const auto s0j = shapeJ(j);
            // assert(s0j == rshapeJ(j));

            for (int k = -1; k <= 1; ++k) {
               const auto s0k = shapeK(k);
               // assert(s0k == rshapeK(k));
               result += s0i * s0j * s0k * By(ci + i, cj + j, ck + k);
            } // end for(k)
         } // end for(j)
      } // end for(i)

      return result;
   }

   auto BzInterp(const auto& Bz, const auto& p0, const auto& cids) {
      using jitTSC = interp::Jit<interp::TSC>;
      using jitCIC = interp::Jit<interp::CIC>;

      const auto [x0, y0, z0] = p0;
      const auto [ci, cj, ck] = cids;

      const jitTSC shapeI(x0 - ci);
      const jitCIC shapeJ(y0 - cj);
      const jitTSC shapeK(z0 - ck);

      auto result = 0.0_fp;
      for (int i = -1; i <= 1; ++i) {
         const auto s0i = shapeI(i);
         for (int j = -1; j <= -1; ++j) {
            const auto s0j = shapeJ(j);
            for (int k = -1; k <= 1; ++k) {
               const auto s0k = shapeK(k);
               result += s0i * s0j * s0k * Bz(ci + i, cj + j, ck + k);
            } // end for(k)
         } // end for(j)
      } // end for(i)

      // assert(result == tmp);
      return result;
   }

   template<int D, bool isB>
   auto FieldToParticleInterp(const auto& F, const auto& p0, const auto& cids) {
      using ShapeFunc = interp::Cached<interp::TSC>;

      static constexpr vec3 b0{-1, -1, -1};
      static constexpr vec3 b1 = interp::rotateOrigin<D>(vec3{1, (D == 1) != isB ? -1 : 0, 1});

      // rotate elements to match loop structure
      const auto [ci, cj, ck] = interp::rotateOrigin<D>(cids);
      const auto [x0, y0, z0] = interp::rotateOrigin<D>(p0);

      const ShapeFunc shapeI(x0);
      const ShapeFunc shapeJ(y0);
      const ShapeFunc shapeK(z0);

      auto result = 0.0_fp;
      for (int i = b0[0]; i <= b1[0]; ++i) {
         const auto s0i = shapeI(i);
         for (int j = b0[1]; j <= b1[1]; ++j) {
            const auto s0j = shapeJ(j);
            for (int k = b0[2]; k <= b1[2]; ++k) {
               const auto s0k = shapeK(k);
               // undo the rotation to get proper indices back
               const auto [x, y, z] = interp::rotateOrigin<D == 2 ? D : !D>(vec3{i + ci, j + cj, k + ck});
               result += s0i * s0j * s0k * F(x, y, z);
            } // end for(k)
         } // end for(j)
      } // end for(i)
      return result;
   } // end FieldToParticle()

   static std::array<double, 6> FieldAtParticle(Particle& p, const auto& emdata) {
      const auto cids = getCIDs(p.location + 0.5);
      const auto p0 = p.location - cids.as_type<compute_t>();

      const auto exc = FieldToParticleInterp<0, 0>(emdata.Ex, p0, cids);
      const auto eyc = FieldToParticleInterp<1, 0>(emdata.Ey, p0, cids);
      const auto ezc = FieldToParticleInterp<2, 0>(emdata.Ez, p0, cids);
      const auto bxc = FieldToParticleInterp<0, 1>(emdata.Bx, p0, cids);
      const auto byc = FieldToParticleInterp<1, 1>(emdata.By, p0, cids);
      const auto bzc = FieldToParticleInterp<2, 1>(emdata.Bz, p0, cids);

      return {exc, eyc, ezc, bxc, byc, bzc};
   } // end FieldAtParticle

   struct BorisPush {
      using emdata_t = electromagnetics::EMData;
      using group_t = ParticleGroup;

      static constexpr std::size_t BC_DEPTH = PMLDepth / 2;

      static void update_velocity(Particle& p, const emdata_t& emdata, auto qdt) {
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

      // #pragma omp declare simd notinbranch
      static void update_position(Particle& p) {
         static constexpr vec3 delta_inv{dt / dx, dt / dy, dt / dz};
         if (p.disabled) { return; }

         const auto new_loc = p.location + delta_inv * p.velocity / p.gamma;
         const auto [inew, jnew, knew] = getCIDs(new_loc);

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
         // Easiest way is to just copy the velocity update and make qdt negative
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
