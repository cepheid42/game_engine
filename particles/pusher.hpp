#ifndef PUSHER_HPP
#define PUSHER_HPP

#include "program_params.hpp"
#include "particles.hpp"
#include "em_data.hpp"
#include "constants.hpp"
#include "interpolation.hpp"

#include <cmath>


namespace tf::particles {
template <int D, typename Strategy>
auto FieldToParticleInterp(const auto& F,
                           const auto& shapeI, const auto& shapeJ, const auto& shapeK,
                           const auto& ci, const auto& cj, const auto& ck)
{
   using IShape = typename Strategy::OuterShape;
   using JShape = typename Strategy::MiddleShape;
   using KShape = typename Strategy::InnerShape;

   auto result = 0.0;
   // auto sum = 0.0;
   for (int i = IShape::Begin; i <= IShape::End; ++i) {
      const auto& s0i = shapeI[i - IShape::Begin];
      for (int j = JShape::Begin; j <= JShape::End; ++j) {
         const auto& s0j = shapeJ[j - JShape::Begin];
         for (int k = KShape::Begin; k <= KShape::End; ++k) {
            const auto& s0k = shapeK[k - KShape::Begin];
            const auto& [x, y, z] = interp::rotateOrigin<D == 2 ? D : !D>(ci + i, cj + j, ck + k);
            result += s0i * s0j * s0k * F(x, y, z);
            // sum += s0i * s0j * s0k;
            // std::println("{} * {} * {} = {}", s0i, s0j, s0k, sum);
         } // end for(k)
      } // end for(j)
   } // end for(i)

   // if (sum != 1.0) {
   //    std::println("Sum: {}", sum);
   //    exit(0);
   // }
   return result;
} // end FieldToParticle()


static std::array<vec3<double>, 2> fieldAtParticle(const Particle& p, const auto& emdata, const auto qdt) {
   using AssShape = interp::InterpolationShape<interpolation_order>;
   using RedShape = interp::InterpolationShape<interpolation_order - 1>;
   using EStrategy = interp::InterpolationStrategy<AssShape, AssShape, RedShape>;
   using BStrategy = interp::InterpolationStrategy<RedShape, RedShape, AssShape>;

   // constexpr auto offset = interpolation_order == 2 ? 0.5 : 1.0;
   // const vec3 loc_full = getCellIndices<double>(p.location + offset);
   // const vec3 loc_half = getCellIndices<double>(p.location + offset + 0.5) - 0.5;

   // First Order
   const vec3 loc_full = getCellIndices<double>(p.location + 1.0);
   const vec3 loc_half = getCellIndices<double>(p.location + 1.0) - 0.5;

   // // Second Order
   // const vec3 loc_full = getCellIndices<double>(p.location + 0.5);
   // const vec3 loc_half = getCellIndices<double>(p.location + 1.5) - 0.5;


   const vec3 fid = loc_full.as_type<std::size_t>();
   const vec3 hid = loc_half.as_type<std::size_t>();

   const vec3 p_full = p.location - loc_full;
   const vec3 p_half = p.location - loc_half;

   const auto xh_r = RedShape::Type::shape_array(p_half[0]);
   const auto yh_r = RedShape::Type::shape_array(p_half[1]);
   const auto zh_r = RedShape::Type::shape_array(p_half[2]);

   const auto xf_a = AssShape::Type::shape_array(p_full[0]);
   const auto yf_a = AssShape::Type::shape_array(p_full[1]);
   const auto zf_a = AssShape::Type::shape_array(p_full[2]);

   // std::println("{}", p.location);
   // std::println("Lfull: {}, Lhalf: {}", loc_full, loc_half);
   // std::println("Pfull: {}, Phalf: {}", p_full, p_half);

   // // First Order
   // // std::println("xh_r: {}", xh_r[0]);
   // std::println("yh_r: {}", yh_r[0]);
   // // std::println("zh_r: {}", zh_r[0]);
   // std::println("xf_a: {}, {}", xf_a[0], xf_a[1]);
   // // std::println("yf_a: {}, {}", yf_a[0], yf_a[1]);
   // std::println("zf_a: {}, {}", zf_a[0], zf_a[1]);

   // // Second Order
   // std::println("xh_r: {}, {}", xh_r[0], xh_r[1]);
   // // std::println("yh_r: {}, {}", yh_r[0], yh_r[1]);
   // std::println("zh_r: {}, {}", zh_r[0], zh_r[1]);
   // // std::println("xf_a: {}, {}, {}", xf_a[0], xf_a[1], xf_a[2]);
   // std::println("yf_a: {}, {}, {}", yf_a[0], yf_a[1], yf_a[2]);
   // // std::println("zf_a: {}, {}, {}", zf_a[0], zf_a[1], zf_a[2]);

   const auto exc = FieldToParticleInterp<0, EStrategy>(emdata.Ex_total, yf_a, zf_a, xh_r, fid[1], fid[2], hid[0]);
   const auto eyc = FieldToParticleInterp<1, EStrategy>(emdata.Ey_total, zf_a, xf_a, yh_r, fid[2], fid[0], hid[1]);
   const auto ezc = FieldToParticleInterp<2, EStrategy>(emdata.Ez_total, xf_a, yf_a, zh_r, fid[0], fid[1], hid[2]);
   const auto bxc = FieldToParticleInterp<0, BStrategy>(emdata.Bx_total, yh_r, zh_r, xf_a, hid[1], hid[2], fid[0]);
   const auto byc = FieldToParticleInterp<1, BStrategy>(emdata.By_total, zh_r, xh_r, yf_a, hid[2], hid[0], fid[1]);
   const auto bzc = FieldToParticleInterp<2, BStrategy>(emdata.Bz_total, xh_r, yh_r, zf_a, hid[0], hid[1], fid[2]);
   // return {};
   return {qdt * vec3{exc, eyc, ezc}, qdt * vec3{bxc, byc, bzc}};
} // end FieldAtParticle

struct BorisPush {
   using emdata_t = electromagnetics::EMData;
   using group_t = ParticleGroup;

   static constexpr std::size_t BC_DEPTH = 3zu;

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

      // Periodic particle BCs
      static constexpr auto fnx = static_cast<double>(Nx - 1);
      static constexpr auto fny = static_cast<double>(Ny - 1);
      static constexpr auto fnz = static_cast<double>(Nz - 1);
      auto new_loc = p.location + delta_inv * p.velocity;
      auto old_loc = p.location;
      if (new_loc[0] <= NHalo or new_loc[0] >= fnx - NHalo) {
         new_loc[0] = fnx + new_loc[0] - 2.0 * std::floor(new_loc[0] + 0.5);
         old_loc[0] = fnx + old_loc[0] - 2.0 * std::floor(old_loc[0] + 0.5);
      }

      if (new_loc[1] <= NHalo or new_loc[1] >= fny - NHalo) {
         new_loc[1] = fny + new_loc[1] - 2.0 * std::floor(new_loc[1] + 0.5);
         old_loc[1] = fny + old_loc[1] - 2.0 * std::floor(old_loc[1] + 0.5);
      }

      if (new_loc[2] <= NHalo or new_loc[2] >= fnz - NHalo) {
         new_loc[2] = fnz + new_loc[2] - 2.0 * std::floor(new_loc[2] + 0.5);
         old_loc[2] = fnz + old_loc[2] - 2.0 * std::floor(old_loc[2] + 0.5);
      }

      // // Outflow particle BCs
      // const auto new_loc = p.location + delta_inv * p.velocity;
      // const auto old_loc = p.location;
      // const auto& [inew, jnew, knew] = getCellIndices(new_loc);
      // p.disabled = inew < BC_DEPTH or inew > Nx - 1 - BC_DEPTH or
      //              jnew < BC_DEPTH or jnew > Ny - 1 - BC_DEPTH or
      //              knew < BC_DEPTH or knew > Nz - 1 - BC_DEPTH;

      p.old_location = old_loc;
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
      // advance_velocity(g, emdata);
      advance_position(g);
      if (step % group_t::SORT_INTERVAL == 0) {
         g.sort_particles();
      }
   } // end operator()
}; // end struct BorisPush
} // end namespace tf::particles


#endif //PUSHER_HPP
