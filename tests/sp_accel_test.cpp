#include "program_params.hpp"
#include "constants.hpp"
#include "em_solver.hpp"
#include "particles.hpp"
#include "pusher.hpp"
#include "vec3.hpp"

#include <cassert>
#include <cmath>
#include <vector>
#include <algorithm>
#include <print>
#include <numeric>

using namespace tf;
using namespace tf::electromagnetics;
using namespace tf::particles;

inline constexpr auto m_e = constants::m_e<compute_t>;
inline constexpr auto q_e = constants::q_e<compute_t>;
inline constexpr auto m_p = constants::m_p<compute_t>;

void location2File(const auto& x, const auto& z, const std::string& path) {
   std::ofstream file(path);
   file << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10);
   for (int i = 0; i < Nt; i++) {
      file << x[i] << ", " << 0.0 << ", " << z[i] << std::endl;
   }
   file.close();
}


double computeMAE(const auto& x, const auto& mu) {
   auto error = std::transform_reduce(x.begin(), x.end(), mu.begin(),
                0.0, std::plus{}, [](const auto a, const auto b) { return std::abs(a - b); });
   return error / static_cast<compute_t>(Nt);
}

template<int D>
double test_pusher(ParticleGroup& g1, const auto& target, auto&& apply_fields) {
   static constexpr std::array min_range = {x_range[0], y_range[0], z_range[0]};
   static constexpr std::array deltas = {dx, dy, dz};

   EMSolver emsolver(Nx, Ny, Nz, cfl, dt);
   apply_fields(emsolver.emdata);
   emsolver.particle_correction();
   // BorisPush::backstep_velocity(g1, emsolver.emdata);
   std::vector<compute_t> pos(Nt);
   pos[0] = min_range[D] + deltas[D] * g1.particles[0].location[D];
   for (std::size_t n = 1; n < Nt; n++) {
      // const auto& p1 = g1.particles[0];
      // std::println("{} {}", p1.velocity, p1.location);
      // no need to advance fields, since they don't change in time
      BorisPush::advance(g1, emsolver.emdata, 1zu);
      pos[n] = min_range[D] + deltas[D] * g1.particles[0].location[D];

      // if (n == 5) { exit(0); }
   }
   return computeMAE(pos, target);
}

void test_ex(const auto E_field, auto&& apply_fields) {
   ParticleGroup g1("electrons", m_e, -q_e, 0);
   g1.particles.emplace_back(
      vec3{2.5_fp, 0.5_fp, 5.5_fp}, // location
      vec3{2.5_fp, 0.5_fp, 5.5_fp}, // old location
      vec3{0.0_fp, 0.0_fp, 0.0_fp}, // velocity
      1.0_fp,                       // weight
      1.0_fp,                       // gamma
      false                         // disabled
   );
   g1.initial_y_position = 0.5_fp;
   const auto ts = math::linspace(0.0, t_end, Nt);
   std::vector<compute_t> xs(Nt);
   std::ranges::transform(ts, xs.begin(), [&](const auto& t) { return 0.5_fp * (-q_e * E_field / m_e) * t * t; });
   const auto error = test_pusher<0>(g1, xs, std::forward<decltype(apply_fields)>(apply_fields));
   std::println("Ex Acceleration Test:\n\tMAE = {:.3e} (≤ 2.202e-11)", error);
   assert(error <= 5.863e-11);
}

void test_ez(const auto E_field, auto&& apply_fields) {
   ParticleGroup g1("electrons", m_e, -q_e, 0);
   g1.particles.emplace_back(
      vec3{5.5_fp, 0.5_fp, 2.5_fp}, // location
      vec3{5.5_fp, 0.5_fp, 2.5_fp}, // old location
      vec3{0.0_fp, 0.0_fp, 0.0_fp}, // velocity
      1.0_fp,                       // weight
      1.0_fp,                       // gamma
      false                         // disabled
   );
   g1.initial_y_position = 0.5_fp;
   const auto ts = math::linspace(0.0, t_end, Nt);
   std::vector<compute_t> xs(Nt);
   std::ranges::transform(ts, xs.begin(), [&](const auto& t) { return 0.5_fp * (-q_e * E_field / m_e) * t * t; });
   const auto error = test_pusher<2>(g1, xs, std::forward<decltype(apply_fields)>(apply_fields));
   std::println("Ez Acceleration Test:\n\tMAE = {:.3e} (≤ 2.202e-11)", error);
   assert(error <= 5.863e-11);
}

void test_ex_quadratic() {
   constexpr auto v_init = 100.0;
   ParticleGroup g1("electrons", m_e, -q_e, 0);
   g1.particles.emplace_back(
      vec3{5.5_fp, 0.5_fp, 25.5_fp}, // location
      vec3{5.5_fp, 0.5_fp, 25.5_fp}, // old location
      vec3{0.0_fp, 0.0_fp, v_init}, // velocity
      1.0_fp,                       // weight
      1.0_fp,                       // gamma
      false                         // disabled
   );
   g1.initial_y_position = 0.5_fp;
   EMSolver emsolver(Nx, Ny, Nz);

   const auto xr = math::linspace(-1.0, 1.0, Nx);
   const auto zr = math::linspace(-1.0, 1.0, Nz);

   for (int i = 0; i < Nx - 1; i++) {
      for (int j = 0; j < Ny; j++) {
         for (int k = 0; k < Nz - 1; k++) {
            const auto r = std::sqrt(math::SQR(xr[i]) + math::SQR(zr[k]));
            const auto theta = std::atan2(zr[k], xr[i]);
            emsolver.emdata.Ex_app(i, j, k) = r * std::cos(theta);
            emsolver.emdata.Ez_app(i, j, k) = r * std::sin(theta);
         }
      }
   }

   emsolver.particle_correction();
   // BorisPush::backstep_velocity(g1, emsolver.emdata);
   std::vector<compute_t> xp(Nt);
   std::vector<compute_t> zp(Nt);
   xp[0] = x_range[0] + dx * g1.particles[0].location[0];
   zp[0] = z_range[0] + dz * g1.particles[0].location[2];
   for (std::size_t n = 1; n < Nt; n++) {
      // const auto& p1 = g1.particles[0];
      // std::println("{} {}", p1.velocity, p1.location);
      // no need to advance fields, since they don't change in time
      BorisPush::advance(g1, emsolver.emdata, 1zu);
      xp[n] = x_range[0] + dx * g1.particles[0].location[0];
      zp[n] = z_range[0] + dz * g1.particles[0].location[2];
   }
   // location2File(xp, zp, "/home/cepheid/TriForce/game_engine/data/larmor.csv");
}

int main() {
   std::println("Single Particle Acceleration Tests.");
   test_ex(-1.0, [](auto& emdata) { emdata.Ex_app.fill(-1.0); });
   test_ez(-1.0, [](auto& emdata) { emdata.Ez_app.fill(-1.0); });

   // test_ex_quadratic();
   return 0;
}