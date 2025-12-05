#ifndef PUSHER_HPP
#define PUSHER_HPP

#include "constants.hpp"
#include "em_data.hpp"
#include "interpolation.hpp"
#include "program_params.hpp"
#include "particles.hpp"
#include "vec3.hpp"

#include <cmath>

namespace tf::particles {
template <int D, typename Strategy>
auto FieldToParticleInterp(const auto& F,
                           const auto& shapeI, const auto& shapeJ, const auto& shapeK,
                           const auto ci, const auto cj, const auto ck)
-> double
{
   using IShape = typename Strategy::OuterShape;
   using JShape = typename Strategy::MiddleShape;
   using KShape = typename Strategy::InnerShape;
   auto result = 0.0;
   for (int i = IShape::Begin; i <= IShape::End; ++i) {
      const auto& s0i = shapeI[i - IShape::Begin];
      for (int j = JShape::Begin; j <= JShape::End; ++j) {
         const auto& s0j = shapeJ[j - JShape::Begin];
         for (int k = KShape::Begin; k <= KShape::End; ++k) {
            const auto& s0k = shapeK[k - KShape::Begin];
            const auto& [x, y, z] = interp::rotateOrigin<D == 2 ? D : !D>(ci + i, cj + j, ck + k);
            result += s0i * s0j * s0k * F[x, y, z];
         } // end for(k)
      } // end for(j)
   } // end for(i)

   // if (result < -2713600.0) {
   //    std::println();
   //    for (int i = IShape::Begin; i <= IShape::End; ++i) {
   //       const auto& s0i = shapeI[i - IShape::Begin];
   //       for (int j = JShape::Begin; j <= JShape::End; ++j) {
   //          const auto& s0j = shapeJ[j - JShape::Begin];
   //          for (int k = KShape::Begin; k <= KShape::End; ++k) {
   //             const auto& s0k = shapeK[k - KShape::Begin];
   //             const auto& [x, y, z] = interp::rotateOrigin<D == 2 ? D : !D>(ci + i, cj + j, ck + k);
   //             std::println("{}, {}, {}: {}, {}, {} * {}", x, y, z, s0i, s0j, s0k, F[x, y, z]);
   //          } // end for(k)
   //       } // end for(j)
   //    } // end for(i)
   //    exit(0);
   // }

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
      XFullShape::Order % 2 == 0 ? 0.5 : 1.0,
      YFullShape::Order % 2 == 0 ? 0.5 : 1.0,
      ZFullShape::Order % 2 == 0 ? 0.5 : 1.0
   };
   const vec3 loc_full = getCellIndices<double>(p.location + offset);
   const vec3 loc_half = loc_full + 0.5;

   const auto& [fx, fy, fz] = loc_full.as_type<std::size_t>();
   const auto& [hx, hy, hz] = loc_half.as_type<std::size_t>();

   const vec3 p_full = p.location - loc_full;
   const vec3 p_half = p.location - loc_half;

   const auto xr =  XRedShape::shape_array(p_half.x);
   const auto yr =  YRedShape::shape_array(p_half.y);
   const auto zr =  ZRedShape::shape_array(p_half.z);
   const auto xa = XFullShape::shape_array(p_full.x);
   const auto ya = YFullShape::shape_array(p_full.y);
   const auto za = ZFullShape::shape_array(p_full.z);

   // todo: add total fields back in eventually
   // const auto exc = FieldToParticleInterp<0, ExStrategy>(emdata.Exf, ya, za, xr, fy, fz, hx);
   const auto eyc = FieldToParticleInterp<1, EyStrategy>(emdata.Eyf, za, xa, yr, fz, fx, hy);
   // const auto ezc = FieldToParticleInterp<2, EzStrategy>(emdata.Ezf, xa, ya, zr, fx, fy, hz);
   // const auto bxc = FieldToParticleInterp<0, BxStrategy>(emdata.Bxf, yr, zr, xa, hy, hz, fx);
   // const auto byc = FieldToParticleInterp<1, ByStrategy>(emdata.Byf, zr, xr, ya, hz, hx, fy);
   // const auto bzc = FieldToParticleInterp<2, BzStrategy>(emdata.Bzf, xr, yr, za, hx, hy, fz);

   // if (eyc * qdt < -4320000.0)
   // {
   //    std::println();
   //    std::println("xa = {}, {}", xa[0], xa[1]);
   //    std::println("ya = {}, {}", ya[0], ya[1]);
   //    std::println("za = {}, {}", za[0], za[1]);
   //    std::println("loc = {}, {}, {}", p.location.x, p.location.y, p.location.z);
   //    std::println("full = {}, {}, {}", loc_full.x, loc_full.y, loc_full.z);
   //    std::println("half = {}, {}, {}", loc_half.x, loc_half.y, loc_half.z);
   //    std::println("p_full = {}, {}, {}", p_full.x, p_full.y, p_full.z);
   //    std::println("p_half = {}, {}, {}", p_half.x, p_half.y, p_half.z);
   //    exit(0);
   // }
   return {qdt * vec3{0.0, eyc, 0.0}, qdt * vec3{0.0, 0.0, 0.0}};
   // return {qdt * vec3{exc, eyc, ezc}, qdt * vec3{bxc, byc, bzc}};
} // end FieldAtParticle

template<ParticleBCType BC>
requires (BC == ParticleBCType::Reflecting)
void apply_particle_bcs(auto& p, const auto& new_loc, const auto& old_loc) {
   // todo: implement this
   (void) p;
   (void) new_loc;
   (void) old_loc;
   assert(false);
}

template<ParticleBCType BC>
requires (BC == ParticleBCType::Periodic)
void apply_particle_bcs(auto& p, auto& new_loc, auto& old_loc) {
   constexpr auto fnx = static_cast<double>(Nx - 1);
   constexpr auto fny = static_cast<double>(Ny - 1);
   constexpr auto fnz = static_cast<double>(Nz - 1);

   if (new_loc[0] <= nHalo or new_loc[0] >= fnx - nHalo) {
      new_loc[0] = fnx + new_loc[0] - 2.0 * std::floor(new_loc[0] + 0.5);
      old_loc[0] = fnx + old_loc[0] - 2.0 * std::floor(old_loc[0] + 0.5);
   }
   if (new_loc[1] <= nHalo or new_loc[1] >= fny - nHalo) {
      new_loc[1] = fny + new_loc[1] - 2.0 * std::floor(new_loc[1] + 0.5);
      old_loc[1] = fny + old_loc[1] - 2.0 * std::floor(old_loc[1] + 0.5);
   }
   if (new_loc[2] <= nHalo or new_loc[2] >= fnz - nHalo) {
      new_loc[2] = fnz + new_loc[2] - 2.0 * std::floor(new_loc[2] + 0.5);
      old_loc[2] = fnz + old_loc[2] - 2.0 * std::floor(old_loc[2] + 0.5);
   }
   p.old_location = old_loc;
   p.location = new_loc;
}

template<ParticleBCType BC>
requires (BC == ParticleBCType::Outflow)
void apply_particle_bcs(auto& p, const auto& new_loc, const auto& old_loc) {
   constexpr std::size_t BC_DEPTH = 3lu;

   const auto& [inew, jnew, knew] = getCellIndices(new_loc);
   const auto disabled = ((inew < BC_DEPTH or inew > Nx - 1 - BC_DEPTH) and !x_collapsed) or
                         // ((jnew < BC_DEPTH or jnew > Ny - 1 - BC_DEPTH) and !y_collapsed) or
                         ((knew < BC_DEPTH or knew > Nz - 1 - BC_DEPTH) and !z_collapsed);
   if (disabled) {
      p.weight = -1.0;
   }
   p.old_location = old_loc;
   p.location = new_loc;
}

struct BorisPush {
   using emdata_t = electromagnetics::emdata_t;
   using group_t = ParticleGroup;

   static void update_velocity(Particle& p, const emdata_t& emdata, const auto qdt) {
      if (p.is_disabled()) { return; }

      const auto& [eps, bet] = fieldAtParticle(p, emdata, qdt);

      // u = gamma * v
      const auto um = p.gamma * p.velocity + eps;
      const auto t = bet / std::sqrt(1.0 + um.length_squared() * constants::over_c_sqr<double>);
      const auto s = 2.0 * t / (1.0 + t.length_squared());
      const auto u_prime = um + cross(um, t);
      const auto u_plus = um + cross(u_prime, s);
      const auto u = u_plus + eps;

      if (std::isnan(calculateGammaV(u)) or std::isnan(u.x) or std::isnan(u.y) or std::isnan(u.z)) {
         std::println();
         std::println("eps: {} | bet: {}", eps, bet);
         std::println("u: {}, {}, {}", u.x, u.y, u.z);
         std::println("vel: {}, {}, {}", p.velocity.x, p.velocity.y, p.velocity.z);
         std::println("loc: {}, {}, {}", p.location.x, p.location.y, p.location.z);
         std::println("old: {}, {}, {}", p.old_location.x, p.old_location.y, p.old_location.z);
         std::println("gamma: {} ({})", p.gamma, calculateGammaV(u));
         exit(0);
      }

      p.gamma = calculateGammaV(u);
      p.velocity = u / p.gamma;
   } // end update_velocity()

   static void update_position(Particle& p) {
      static constexpr vec3 delta_inv{dt / dx, dt / dy, dt / dz};

      if (p.is_disabled()) { return; }

      auto new_loc = p.location + (delta_inv * p.velocity);
      auto old_loc = p.location;

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

   static void advance(auto& g, const auto& emdata) requires(push_enabled) {
      advance_velocity(g, emdata);
      advance_position(g);
      g.cell_map_updated = false;
      g.is_sorted = false;
   } // end advance()

   static void advance(auto&, const auto&) requires (!push_enabled) {}
   static void backstep_velocity(auto&, const auto&) requires (!push_enabled) {}
}; // end struct BorisPush
} // end namespace tf::particles


#endif //PUSHER_HPP
