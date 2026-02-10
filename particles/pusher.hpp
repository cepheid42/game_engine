#ifndef PUSHER_HPPx0
#define PUSHER_HPP

#include "constants.hpp"
#include "em_data.hpp"
#include "interpolation.hpp"
#include "particles.hpp"
#include "program_params.hpp"
#include "vec3.hpp"

#include <cmath>

namespace tf::particles {
template <int D, typename Strategy>
auto FieldToParticleInterp(const auto& F,
                           const auto& shapeI, const auto& shapeJ, const auto& shapeK,
                           const auto ci, const auto cj, const auto ck)
-> double
{
   using IShape = Strategy::OuterShape;
   using JShape = Strategy::MiddleShape;
   using KShape = Strategy::InnerShape;
   auto result = 0.0;
   for (int i = IShape::Begin; i <= IShape::End; ++i) {
      const auto& s0i = shapeI[i - IShape::Begin];
      for (int j = JShape::Begin; j <= JShape::End; ++j) {
         const auto& s0j = shapeJ[j - JShape::Begin];
         for (int k = KShape::Begin; k <= KShape::End; ++k) {
            const auto& s0k = shapeK[k - KShape::Begin];
            const auto& [x, y, z] = interp::rotateOrigin<D == 2 ? D : !D>(ci + i, cj + j, ck + k);
            result += s0i * s0j * s0k * F(x, y, z);
         } // end for(k)
      } // end for(j)
   } // end for(i)
   return result;
} // end FieldToParticle()

static auto fieldAtParticle(const Particle& p, const auto& emdata, const auto qdt)
-> std::array<vec3<double>, 2>
{
   using XFullShape = interp::InterpolationShape<x_collapsed ? 1 : interpolation_order>::Type;
   using YFullShape = interp::InterpolationShape<y_collapsed ? 1 : interpolation_order>::Type;
   using ZFullShape = interp::InterpolationShape<z_collapsed ? 1 : interpolation_order>::Type;
   using XRedShape  = interp::InterpolationShape<x_collapsed ? 0 : interpolation_order - 1>::Type;
   using YRedShape  = interp::InterpolationShape<y_collapsed ? 0 : interpolation_order - 1>::Type;
   using ZRedShape  = interp::InterpolationShape<z_collapsed ? 0 : interpolation_order - 1>::Type;

   using ExStrategy = interp::InterpolationStrategy<YFullShape, ZFullShape, XRedShape>;
   using EyStrategy = interp::InterpolationStrategy<ZFullShape, XFullShape, YRedShape>;
   using EzStrategy = interp::InterpolationStrategy<XFullShape, YFullShape, ZRedShape>;
   using BxStrategy = interp::InterpolationStrategy<YRedShape,  ZRedShape,  XFullShape>;
   using ByStrategy = interp::InterpolationStrategy<ZRedShape,  XRedShape,  YFullShape>;
   using BzStrategy = interp::InterpolationStrategy<XRedShape,  YRedShape,  ZFullShape>;

   static constexpr vec3 offset{
      XFullShape::Order % 2 == 0 ? 0.5f : 1.0f,
      YFullShape::Order % 2 == 0 ? 0.5f : 1.0f,
      ZFullShape::Order % 2 == 0 ? 0.5f : 1.0f
   };
   const vec3 loc_full = getCellIndices<float>(p.location + offset);
   const vec3 loc_half = loc_full + 0.5f;

   const vec3 fid = loc_full.to_uint();
   const vec3 hid = loc_half.to_uint();

   const vec3 p_full = p.location - loc_full;
   const vec3 p_half = p.location - loc_half;

   const auto xh_r =  XRedShape::shape_array(p_half.x);
   const auto yh_r =  YRedShape::shape_array(p_half.y);
   const auto zh_r =  ZRedShape::shape_array(p_half.z);
   const auto xf_a = XFullShape::shape_array(p_full.x);
   const auto yf_a = YFullShape::shape_array(p_full.y);
   const auto zf_a = ZFullShape::shape_array(p_full.z);

   const auto exc = FieldToParticleInterp<0, ExStrategy>(emdata.Ex_total, yf_a, zf_a, xh_r, fid.y, fid.z, hid.x);
   const auto eyc = FieldToParticleInterp<1, EyStrategy>(emdata.Ey_total, zf_a, xf_a, yh_r, fid.z, fid.x, hid.y);
   const auto ezc = FieldToParticleInterp<2, EzStrategy>(emdata.Ez_total, xf_a, yf_a, zh_r, fid.x, fid.y, hid.z);
   const auto bxc = FieldToParticleInterp<0, BxStrategy>(emdata.Bx_total, yh_r, zh_r, xf_a, hid.y, hid.z, fid.x);
   const auto byc = FieldToParticleInterp<1, ByStrategy>(emdata.By_total, zh_r, xh_r, yf_a, hid.z, hid.x, fid.y);
   const auto bzc = FieldToParticleInterp<2, BzStrategy>(emdata.Bz_total, xh_r, yh_r, zf_a, hid.x, hid.y, fid.z);

   return {qdt * vec3{exc, eyc, ezc}, qdt * vec3{bxc, byc, bzc}};
} // end FieldAtParticle

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
void apply_particle_bcs(auto& p, auto& new_loc, auto& old_loc) {
   constexpr auto fnx = static_cast<float>(Nx - 1);
   constexpr auto fny = static_cast<float>(Ny - 1);
   constexpr auto fnz = static_cast<float>(Nz - 1);

   if (new_loc[0] <= nHalo or new_loc[0] >= fnx - nHalo) {
      new_loc[0] = fnx + new_loc[0] - 2.0f * std::floor(new_loc[0] + 0.5f);
      old_loc[0] = fnx + old_loc[0] - 2.0f * std::floor(old_loc[0] + 0.5f);
   }
   if (new_loc[1] <= nHalo or new_loc[1] >= fny - nHalo) {
      new_loc[1] = fny + new_loc[1] - 2.0f * std::floor(new_loc[1] + 0.5f);
      old_loc[1] = fny + old_loc[1] - 2.0f * std::floor(old_loc[1] + 0.5f);
   }
   if (new_loc[2] <= nHalo or new_loc[2] >= fnz - nHalo) {
      new_loc[2] = fnz + new_loc[2] - 2.0f * std::floor(new_loc[2] + 0.5f);
      old_loc[2] = fnz + old_loc[2] - 2.0f * std::floor(old_loc[2] + 0.5f);
   }

   p.old_location = old_loc;
   p.location = new_loc;
}

template<ParticleBCType BC>
requires (BC == ParticleBCType::Outflow)
void apply_particle_bcs(auto& p, const auto& new_loc, const auto& old_loc) {
   constexpr std::size_t BC_DEPTH = 3zu;

   const auto [inew, jnew, knew] = getCellIndices(new_loc);
   const auto disabled = ((inew < BC_DEPTH or inew > Nx - 1 - BC_DEPTH) and !x_collapsed) or
                         ((jnew < BC_DEPTH or jnew > Ny - 1 - BC_DEPTH) and !y_collapsed) or
                         ((knew < BC_DEPTH or knew > Nz - 1 - BC_DEPTH) and !z_collapsed);
   if (disabled) {
      p.weight = -1.0f;
   }
   p.old_location = old_loc;
   p.location = new_loc;
}

struct BorisPush {
   using emdata_t = electromagnetics::EMData;
   using group_t = ParticleGroup;

   static void update_velocity(Particle& p, const emdata_t& emdata, const auto qdt) {
      if (p.is_disabled()) { return; }

      const auto& [eps, bet] = fieldAtParticle(p, emdata, qdt);

      // u = gamma * v
      const auto um = p.gamma * p.velocity + eps;
      const auto t = bet / std::sqrt(1.0 + um.length_squared() * constants::over_c_sqr);
      const auto s = 2.0 * t / (1.0 + t.length_squared());
      const auto u_prime = um + cross(um, t);
      const auto u_plus = um + cross(u_prime, s);
      const auto u = u_plus + eps;

      p.gamma = std::sqrt(1.0 + u.length_squared() * constants::over_c_sqr);
      p.velocity = u / p.gamma;
   } // end update_velocity()

   static void update_position(Particle& p) {
      static constexpr vec3 delta_inv{dt / dx, dt / dy, dt / dz};

      if (p.is_disabled()) { return; }

      auto old_loc = p.location;
      auto new_loc = p.location + (delta_inv * p.velocity).to_float();

      apply_particle_bcs<PBCSelect>(p, new_loc, old_loc);
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

   static void backstep_velocity(auto& g, const auto& emdata) requires(push_enabled) {
      #pragma omp parallel for num_threads(nThreads)
      for (std::size_t pid = 0; pid < g.num_particles(); pid++) {
         update_velocity(g.particles[pid], emdata, -0.5 * g.qdt_over_2m);
      }
   }

   static void advance(auto& g, const auto& emdata, const auto) requires(push_enabled) {
      if (g.is_photons) { return; }
      advance_velocity(g, emdata);
      advance_position(g);
      g.cell_map_updated = false;
      g.is_sorted = false;
   } // end advance()

   static void advance(auto&, const auto&, const auto) requires (!push_enabled) {}
   static void backstep_velocity(auto&, const auto&) requires (!push_enabled) {}
}; // end struct BorisPush
} // end namespace tf::particles


#endif //PUSHER_HPP
