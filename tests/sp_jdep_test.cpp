#include "program_params.hpp"
#include "constants.hpp"
#include "em_solver.hpp"
#include "pusher.hpp"
#include "current_deposition.hpp"
#include "vec3.hpp"

#include <cassert>
#include <cmath>
#include <vector>
#include <algorithm>
#include <print>

using namespace tf;
using namespace tf::electromagnetics;
using namespace tf::particles;

void electron_test_x() {
   constexpr auto m_e = constants::m_e<compute_t>;
   constexpr auto q_e = constants::q_e<compute_t>;
   ParticleGroup g1("electrons", m_e, -q_e, 0);
   constexpr auto init_vel = 250000.0;
   g1.particles.emplace_back(
      vec3{10.5_fp, 0.5_fp, 25.5},    // location
      vec3{10.5_fp, 0.5_fp, 25.5},    // old location
      vec3{init_vel, 0.0_fp, 0.0_fp}, // velocity
      1.0_fp,                         // weight
      1.0_fp,                         // gamma
      false                           // disabled
   );
   g1.initial_y_position = 0.5_fp;

   EMSolver emsolver(Nx, Ny, Nz, cfl, dt);
   constexpr BorisPush particle_push{};
   constexpr CurrentDeposition current_dep{};
   BorisPush::backstep_velocity(g1, emsolver.emdata);

   constexpr auto j_target = -q_e * init_vel;
   std::vector<compute_t> J(Nt);
   for (std::size_t n = 0; n < Nt; n++) {
      emsolver.advance(dt * static_cast<compute_t>(n));
      particle_push(g1, emsolver.emdata, n);
      current_dep(emsolver.emdata, g1);
      const auto& p = g1.particles[0];
      assert(p.velocity[1] == 0.0);
      assert(p.velocity[2] == 0.0);
      const auto jx_sum = std::ranges::fold_left(emsolver.emdata.Jx.data(), 0.0, std::plus<compute_t>{});
      const auto jy_sum = std::ranges::fold_left(emsolver.emdata.Jy.data(), 0.0, std::plus<compute_t>{});
      const auto jz_sum = std::ranges::fold_left(emsolver.emdata.Jz.data(), 0.0, std::plus<compute_t>{});
      J[n] = std::sqrt(math::SQR(jx_sum) + math::SQR(jy_sum) + math::SQR(jz_sum));
   }
   auto error = 0.0_fp;
   for (std::size_t i = 0; i < Nt; i++) {
      error += std::abs(J[i] - j_target);
   }
   error /= static_cast<compute_t>(Nt);
   std::println("Electron X-Deposition Test: MAE = {:.3e} <= 2.5e-7", error);
   assert(error <= 2.5e-7);
}

void electron_test_z() {
   constexpr auto m_e = constants::m_e<compute_t>;
   constexpr auto q_e = constants::q_e<compute_t>;
   ParticleGroup g1("electrons", m_e, -q_e, 0);
   constexpr auto init_vel = 250000.0;
   g1.particles.emplace_back(
      vec3{25.5_fp, 0.5_fp, 10.5},    // location
      vec3{25.5_fp, 0.5_fp, 10.5},    // old location
      vec3{0.0_fp, 0.0_fp, init_vel}, // velocity
      1.0_fp,                         // weight
      1.0_fp,                         // gamma
      false                           // disabled
   );
   g1.initial_y_position = 0.5_fp;

   EMSolver emsolver(Nx, Ny, Nz, cfl, dt);
   constexpr BorisPush particle_push{};
   constexpr CurrentDeposition current_dep{};
   BorisPush::backstep_velocity(g1, emsolver.emdata);

   constexpr auto j_target = -q_e * init_vel;
   std::vector<compute_t> J(Nt);
   for (std::size_t n = 0; n < Nt; n++) {
      emsolver.advance(dt * static_cast<compute_t>(n));
      particle_push(g1, emsolver.emdata, n);
      current_dep(emsolver.emdata, g1);
      const auto& p = g1.particles[0];
      assert(p.velocity[0] == 0.0);
      assert(p.velocity[1] == 0.0);
      const auto jx_sum = std::ranges::fold_left(emsolver.emdata.Jx.data(), 0.0, std::plus<compute_t>{});
      const auto jy_sum = std::ranges::fold_left(emsolver.emdata.Jy.data(), 0.0, std::plus<compute_t>{});
      const auto jz_sum = std::ranges::fold_left(emsolver.emdata.Jz.data(), 0.0, std::plus<compute_t>{});
      J[n] = std::sqrt(math::SQR(jx_sum) + math::SQR(jy_sum) + math::SQR(jz_sum));
   }
   auto error = 0.0_fp;
   for (std::size_t i = 0; i < Nt; i++) {
      error += std::abs(J[i] - j_target);
   }
   error /= static_cast<compute_t>(Nt);
   std::println("Electron Z-Deposition Test: MAE = {:.3e} <= 2.5e-7", error);
   assert(error <= 2.5e-7);
}



int main() {
   std::println("Single Particle Current Deposition Tests.");
   electron_test_x();
   electron_test_z();
   return 0;
}