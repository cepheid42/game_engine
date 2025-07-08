#include "program_params.hpp"
#include "constants.hpp"
#include "em_solver.hpp"
#include "pusher.hpp"
#include "current_deposition.hpp"
#include "vec3.hpp"
#include "metrics.hpp"

#include <cassert>
#include <vector>
#include <print>

using namespace tf;
using namespace tf::electromagnetics;
using namespace tf::particles;
using namespace tf::metrics;

inline constexpr auto m_e = constants::m_e<compute_t>;
inline constexpr auto q_e = constants::q_e<double> / 1.0e8;
constexpr auto v_init = 3.9945658078880e1;


double test_deposition(ParticleGroup& g1, ParticleGroup& g2) {
   EMSolver emsolver(Nx, Ny, Nz, cfl, dt);
   emsolver.particle_correction();

   BorisPush::backstep_velocity(g1, emsolver.emdata);
   BorisPush::backstep_velocity(g2, emsolver.emdata);

   std::vector<compute_t> J(Nt);
   for (std::size_t n = 0; n < Nt; n++) {
      // const auto& p1 = g1.particles[0];
      // std::println("{} {}", p1.velocity, p1.location);

      emsolver.advance(static_cast<double>(n) * dt);

      BorisPush::advance(g1, emsolver.emdata, 1zu);
      BorisPush::advance(g2, emsolver.emdata, 1zu);

      CurrentDeposition::advance(g1, emsolver.emdata);
      CurrentDeposition::advance(g2, emsolver.emdata);

      const auto jx_sum = std::ranges::fold_left(emsolver.emdata.Jx.vec_data(), 0.0_fp, std::plus<compute_t>{});
      const auto jy_sum = std::ranges::fold_left(emsolver.emdata.Jy.vec_data(), 0.0_fp, std::plus<compute_t>{});
      const auto jz_sum = std::ranges::fold_left(emsolver.emdata.Jz.vec_data(), 0.0_fp, std::plus<compute_t>{});
      J[n] = jx_sum + jy_sum + jz_sum;
   }

   constexpr auto J_target = 0.0;
   auto error = 0.0_fp;
   for (std::size_t i = 0; i < Nt; i++) { error += std::abs(J[i] - J_target); }
   return error / static_cast<compute_t>(Nt);
}

void electron_test_x() {
   ParticleGroup g1("electrons", m_e, -q_e, 0);
   g1.particles.emplace_back(
      vec3{2.5_fp, 0.5_fp, 5.5_fp}, // location
      vec3{2.5_fp, 0.5_fp, 5.5_fp}, // old location
      vec3{v_init, 0.0_fp, 0.0_fp}, // velocity
      1.0_fp,                       // weight
      1.0_fp,                       // gamma
      false                         // disabled
   );
   g1.initial_y_position = 0.5_fp;

   ParticleGroup g2("positrons", m_e, +q_e, 0);
   g2.particles.emplace_back(
      vec3{2.5_fp, 0.5_fp, 5.5_fp}, // location
      vec3{2.5_fp, 0.5_fp, 5.5_fp}, // old location
      vec3{v_init, 0.0_fp, 0.0_fp}, // velocity
      1.0_fp,                       // weight
      1.0_fp,                       // gamma
      false                         // disabled
   );
   g2.initial_y_position = 0.5_fp;

   const auto error = test_deposition(g1, g2);
   std::println("Electron/Positron X-direction JDep\n\tMAE = {:.3e} (≤ 1.0e-16)", error);
   assert(error <= 1.0e-16);
}

void electron_test_z() {
   ParticleGroup g1("electrons", m_e, -q_e, 0);
   g1.particles.emplace_back(
      vec3{5.5_fp, 0.5_fp, 2.5_fp}, // location
      vec3{5.5_fp, 0.5_fp, 2.5_fp}, // old location
      vec3{0.0_fp, 0.0_fp, v_init}, // velocity
      1.0_fp,                       // weight
      1.0_fp,                       // gamma
      false                         // disabled
   );
   g1.initial_y_position = 0.5_fp;

   ParticleGroup g2("positrons", m_e, +q_e, 0);
   g2.particles.emplace_back(
      vec3{5.5_fp, 0.5_fp, 2.5_fp}, // location
      vec3{5.5_fp, 0.5_fp, 2.5_fp}, // old location
      vec3{0.0_fp, 0.0_fp, v_init}, // velocity
      1.0_fp,                       // weight
      1.0_fp,                       // gamma
      false                         // disabled
   );
   g2.initial_y_position = 0.5_fp;

   const auto error = test_deposition(g1, g2);
   std::println("Electron/Positron Z-direction JDep\n\tMAE = {:.3e} (≤ 1.0e-16)", error);
   assert(error <= 1.0e-16);
}

int main() {
   std::println("Single Particle Current Deposition Tests.");
   electron_test_x();
   electron_test_z();
   return 0;
}