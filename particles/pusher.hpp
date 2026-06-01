#ifndef PUSHER_HPP
#define PUSHER_HPP

#include "constants.hpp"
#include "em_data.hpp"
#include "interpolation.hpp"
#include "particles.hpp"
#include "program_params.hpp"
#include "vec3.hpp"

#include <cmath>

namespace tf::particles {
using XFullShape = interp::InterpolationShape<x_collapsed ? 1 : interpolation_order>::Type;
using YFullShape = interp::InterpolationShape<y_collapsed ? 1 : interpolation_order>::Type;
using ZFullShape = interp::InterpolationShape<z_collapsed ? 1 : interpolation_order>::Type;
using XRedShape  = interp::InterpolationShape<x_collapsed ? 0 : interpolation_order - 1>::Type;
using YRedShape  = interp::InterpolationShape<y_collapsed ? 0 : interpolation_order - 1>::Type;
using ZRedShape  = interp::InterpolationShape<z_collapsed ? 0 : interpolation_order - 1>::Type;

using ExStrategy = interp::InterpolationStrategy<0, YFullShape, ZFullShape, XRedShape>;
using EyStrategy = interp::InterpolationStrategy<1, ZFullShape, XFullShape, YRedShape>;
using EzStrategy = interp::InterpolationStrategy<2, XFullShape, YFullShape, ZRedShape>;
using BxStrategy = interp::InterpolationStrategy<0, YRedShape,  ZRedShape,  XFullShape>;
using ByStrategy = interp::InterpolationStrategy<1, ZRedShape,  XRedShape,  YFullShape>;
using BzStrategy = interp::InterpolationStrategy<2, XRedShape,  YRedShape,  ZFullShape>;


template<ParticleBCType BC>
requires (BC == ParticleBCType::Reflecting)
void apply_particle_bcs(auto& p, const auto& new_loc, const auto& old_loc) {
   (void) p;
   (void) new_loc;
   (void) old_loc;
   assert(false);
}

template<ParticleBCType BC>
requires (BC == ParticleBCType::Periodic)
void apply_particle_bcs(auto& p) {
   constexpr auto fnx = static_cast<double>(Nx - 1);
   constexpr auto fny = static_cast<double>(Ny - 1);
   constexpr auto fnz = static_cast<double>(Nz - 1);

   if (p.location[0] <= PBCDepth or p.location[0] >= fnx - PBCDepth) {
      p.location[0] = fnx + p.location[0] - 2.0 * std::floor(p.location[0] + 0.5);
      p.old_location[0] = fnx + p.old_location[0] - 2.0 * std::floor(p.old_location[0] + 0.5);
   }
   if (p.location[1] <= PBCDepth or p.location[1] >= fny - PBCDepth) {
      p.location[1] = fny + p.location[1] - 2.0 * std::floor(p.location[1] + 0.5);
      p.old_location[1] = fny + p.old_location[1] - 2.0 * std::floor(p.old_location[1] + 0.5);
   }
   if (p.location[2] <= PBCDepth or p.location[2] >= fnz - PBCDepth) {
      p.location[2] = fnz + p.location[2] - 2.0 * std::floor(p.location[2] + 0.5);
      p.old_location[2] = fnz + p.old_location[2] - 2.0 * std::floor(p.old_location[2] + 0.5);
   }
}

template<ParticleBCType BC>
requires (BC == ParticleBCType::Outflow)
void apply_particle_bcs(auto& p) {
   const auto [inew, jnew, knew] = getCellIndices(p.location);
   const auto disabled = (!x_collapsed and (inew < PBCDepth or inew > Nx - 2 - PBCDepth)) or
                         (!y_collapsed and (jnew < PBCDepth or jnew > Ny - 2 - PBCDepth)) or
                         (!z_collapsed and (knew < PBCDepth or knew > Nz - 2 - PBCDepth));
   if (disabled) {
      p.weight = -1.0f;
   }
}

template <typename Strategy>
auto FieldToParticleInterp(const auto& F,
                           const auto& shapeI, const auto& shapeJ, const auto& shapeK,
                           const auto ci, const auto cj, const auto ck)
-> double
{
   using IShape = Strategy::OuterShape;
   using JShape = Strategy::MiddleShape;
   using KShape = Strategy::InnerShape;
   auto result = 0.0;
   for (auto i = IShape::Begin; i <= IShape::End; ++i) {
      const auto& s0i = shapeI[i];
      for (auto j = JShape::Begin; j <= JShape::End; ++j) {
         const auto& s0j = shapeJ[j];
         for (auto k = KShape::Begin; k <= KShape::End; ++k) {
            const auto& s0k = shapeK[k];
            const auto& [x, y, z] = Strategy::permute(ci + i, cj + j, ck + k);
            result += s0i * s0j * s0k * F(x, y, z);
         } // end for(k)
      } // end for(j)
   } // end for(i)
   return result;
} // end FieldToParticle()

auto fieldAtParticle(Particle& p, const auto& emdata, const auto qdt)
-> std::array<vec3<double>, 2>
{
   static constexpr vec3 offset{
      interpolation_order == 2 ? 0.5 : 0.0,
      interpolation_order == 2 ? 0.5 : 0.0,
      interpolation_order == 2 ? 0.5 : 0.0
   };

   const vec3 loc_full = getCellIndices<double>(p.location);
   const vec3 loc_half = getCellIndices<double>(p.location - offset);

   const vec3 fid = loc_full.to_uint();
   const vec3 hid = loc_half.to_uint();

   const vec3 p_full = p.location - loc_full;
   const vec3 p_half = p.location - loc_half - 0.5;

   const auto xh =  XRedShape::shape_array(p_half.x);
   const auto yh =  YRedShape::shape_array(p_half.y);
   const auto zh =  ZRedShape::shape_array(p_half.z);
   const auto xf = XFullShape::shape_array(p_full.x);
   const auto yf = YFullShape::shape_array(p_full.y);
   const auto zf = ZFullShape::shape_array(p_full.z);

   const auto exc = FieldToParticleInterp<ExStrategy>(emdata.Ex_total, yf, zf, xh, fid.y, fid.z, hid.x);
   const auto eyc = FieldToParticleInterp<EyStrategy>(emdata.Ey_total, zf, xf, yh, fid.z, fid.x, hid.y);
   const auto ezc = FieldToParticleInterp<EzStrategy>(emdata.Ez_total, xf, yf, zh, fid.x, fid.y, hid.z);

   const auto bxc = FieldToParticleInterp<BxStrategy>(emdata.Bx_total, yh, zh, xf, hid.y, hid.z, fid.x);
   const auto byc = FieldToParticleInterp<ByStrategy>(emdata.By_total, zh, xh, yf, hid.z, hid.x, fid.y);
   const auto bzc = FieldToParticleInterp<BzStrategy>(emdata.Bz_total, xh, yh, zf, hid.x, hid.y, fid.z);

   return {(qdt / constants::c) * vec3{exc, eyc, ezc}, qdt * vec3{bxc, byc, bzc}};
} // end FieldAtParticle

template<ParticlePushType P>
struct ParticleVelocityUpdate {
   static void operator()(Particle&, const auto&, auto)
   requires (P == ParticlePushType::Ballistic)
   {} // end Ballistic Velocity Update

   static void operator()(Particle& p, const auto& emdata, const auto qdt)
   requires (P == ParticlePushType::Boris)
   {
      // todo: maybe figure a way to interpolate all Ex first for every particle, then Ey, etc...
      //       might help improve performance since it's not trying to shuffle all the fields through
      //       at the same time.
      const auto& [eps, bet] = fieldAtParticle(p, emdata, qdt);
      const auto um = p.beta_gamma + eps;
      const auto t = bet / std::sqrt(1.0 + um.length_squared());
      const auto s = 2.0 * t / (1.0 + t.length_squared());
      const auto u_prime = um + cross_product(um, t);
      const auto u_plus = um + cross_product(u_prime, s);
      const auto u = u_plus + eps;
      p.beta_gamma = u;
   } // end Boris Velocity Update

   static void operator()(Particle& p, const auto& emdata, const auto qdt)
   requires (P == ParticlePushType::HigueraCary)
   {
      const auto& [eps, bet] = fieldAtParticle(p, emdata, qdt);

      const auto um = p.beta_gamma + eps;
      const auto tau2 = bet.length_squared();
      const auto gamma_m2 = 1.0 + um.length_squared();
      const auto sigma = gamma_m2 - tau2;
      const auto u_star2 = math::SQR(dot_product(um, bet / constants::c));
      const auto gamma_new = std::sqrt(0.5 * (sigma + std::sqrt(math::SQR(sigma) + 4.0 * (tau2 + u_star2))));
      const auto t = bet / gamma_new;
      const auto s = 1.0 / (1.0 + t.length_squared());
      const auto up = s * (um + dot_product(um, t) * t + cross_product(um, t));
      const auto u = up + eps + cross_product(up, t);
      p.beta_gamma = u;
   } // end Higuera-Cary Velocity Update
}; // end struct ParticleVelocityUpdate

struct ParticlePusher {
   using emdata_t = electromagnetics::EMData;
   using group_t = ParticleGroup;

   static constexpr vec3 cdt_delta_inv{constants::c * dt / dx, constants::c * dt / dy, constants::c * dt / dz};

   static void update_position(Particle& p) {
      p.old_location = p.location;
      p.location += (cdt_delta_inv * p.beta_gamma / p.gamma());
      apply_particle_bcs<PBCSelect>(p);
   } // end update_position()

   static void advance_position(group_t& g) {
      #pragma omp parallel for num_threads(nThreads)
      for (auto pid = 0zu; pid < g.num_particles(); pid++) {
         if (g.particles[pid].is_disabled()) { continue; }
         update_position(g.particles[pid]);
      }
   } // end first_advance_position

   static void advance_velocity(group_t& g, const emdata_t& emdata) {
      #pragma omp parallel for num_threads(nThreads)
      for (auto pid = 0zu; pid < g.num_particles(); pid++) {
         if (g.particles[pid].is_disabled()) { continue; }
         ParticleVelocityUpdate<ParticlePushSelect>()(g.particles[pid], emdata, g.qdt_over_2m);
      }
   } // end advance_velocity

   static void backstep_velocity(group_t& g, const emdata_t& emdata) {
      #pragma omp parallel for num_threads(nThreads)
      for (std::size_t pid = 0; pid < g.num_particles(); pid++) {
         ParticleVelocityUpdate<ParticlePushSelect>()(g.particles[pid], emdata, -0.5 * g.qdt_over_2m);
      }
   }

   static void advance(auto& g, const auto& emdata, const auto step) requires(push_enabled) {
      if (g.is_photons) { return; }
      g.reset_positions();

      if (step % sort_frequency == 0) { g.sort_particles(); }

      advance_velocity(g, emdata); // n -> n+1
      advance_position(g);         // n-1/2 -> n+1/2

      g.cell_map_updated = false;
      g.is_sorted = false;
   } // end advance()

   static void advance(auto&, const auto&, const auto) requires (!push_enabled) {}
   static void backstep_velocity(auto&, const auto&) requires (!push_enabled) {}
}; // end struct BorisPush
} // end namespace tf::particles


#endif //PUSHER_HPP
