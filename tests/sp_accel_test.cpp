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

using namespace tf;
using namespace tf::electromagnetics;
using namespace tf::particles;

void electron_test() {
   constexpr auto m_e = constants::m_e<compute_t>;
   constexpr auto q_e = constants::q_e<compute_t>;
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

   EMSolver emsolver(Nx, Ny, Nz, cfl, dt);
   constexpr auto E_field = -1.0_fp;
   emsolver.emdata.Ex_app.fill(E_field);
   emsolver.particle_correction();

   constexpr BorisPush particle_push{};
   BorisPush::backstep_velocity(g1, emsolver.emdata);

   std::vector<compute_t> pos(Nt);

   for (std::size_t n = 0; n < Nt; n++) {
      particle_push(g1, emsolver.emdata, 1zu);
      g1.reset_y_positions();

      pos[n] = x_range[0] + dx * g1.particles[0].location[0];
   }

   const auto ts = math::linspace(0.0, t_end, Nt);
   std::vector<compute_t> xs(Nt);
   std::ranges::transform(ts, xs.begin(), [&](const auto& t) { return 0.5_fp * (-q_e * E_field / m_e) * t * t; });

   auto error = 0.0_fp;
   for (int i = 0; i < Nt; i++) {
      error += std::abs(xs[i] - pos[i]);
   }
   error /= static_cast<compute_t>(Nt);

   std::println("Electron Acceleration Test: MAE = {:.3e} <= 3.0e-11", error);
   assert(error <= 3.0e-11);
}

void ion_test() {
   constexpr auto m_p = constants::m_p<compute_t>;
   constexpr auto q_p = constants::q_e<compute_t>;
   ParticleGroup g1("ions", m_p, q_p, 1);
   g1.particles.emplace_back(
      vec3{2.5_fp, 0.5_fp, 5.5_fp}, // location
      vec3{2.5_fp, 0.5_fp, 5.5_fp}, // old location
      vec3{0.0_fp, 0.0_fp, 0.0_fp}, // velocity
      1.0_fp,                       // weight
      1.0_fp,                       // gamma
      false                         // disabled
   );

   g1.initial_y_position = 0.5_fp;

   EMSolver emsolver(Nx, Ny, Nz, cfl, dt);
   constexpr auto E_field = 1.0_fp;
   emsolver.emdata.Ex_app.fill(E_field);
   emsolver.particle_correction();

   constexpr BorisPush particle_push{};
   BorisPush::backstep_velocity(g1, emsolver.emdata);

   std::vector<compute_t> pos(Nt);

   for (std::size_t n = 0; n < Nt; n++) {
      particle_push(g1, emsolver.emdata, 1zu);
      g1.reset_y_positions();

      pos[n] = x_range[0] + dx * g1.particles[0].location[0];
   }

   const auto ts = math::linspace(0.0, t_end, Nt);
   std::vector<compute_t> xs(Nt);
   std::ranges::transform(ts, xs.begin(), [&](const auto& t) { return 0.5_fp * (q_p * E_field / m_p) * t * t; });

   auto error = 0.0_fp;
   for (int i = 0; i < Nt; i++) {
      error += std::abs(xs[i] - pos[i]);
   }
   error /= static_cast<compute_t>(Nt);
   std::println("Ion Acceleration Test: MAE = {:.3e} <= 2.0e-14", error);
   assert(error <= 2.0e-14);
}

int main() {
   std::println("Single Particle Acceleration Tests.");
   electron_test();
   ion_test();
   return 0;
}