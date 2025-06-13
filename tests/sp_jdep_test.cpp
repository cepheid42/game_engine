#include "program_params.hpp"
#include "constants.hpp"
#include "em_solver.hpp"
#include "pusher.hpp"
#include "current_deposition.hpp"
#include "vec3.hpp"
#include "metrics.hpp"

#include <cassert>
// #include <cmath>
#include <vector>
// #include <algorithm>
#include <print>

using namespace tf;
using namespace tf::electromagnetics;
using namespace tf::particles;
using namespace tf::metrics;

inline constexpr auto m_e = constants::m_e<compute_t>;
inline constexpr auto q_e = 2.5600008901198904e-13_fp;
inline constexpr auto v_init = 250000.0_fp;
inline constexpr auto J_target = -q_e * v_init / (dx * dy * dz);


Metrics create_metrics(const std::string& dir, EMSolver& em, const ParticleGroup& g1) {
   Metrics metrics(dir);

   metrics.addMetric(
      std::make_unique<EMFieldsMetric>(
         std::unordered_map<std::string, Array3D<compute_t>*>{
            {"Jx", &em.emdata.Jx},
            {"Jy", &em.emdata.Jy},
            {"Jz", &em.emdata.Jz}
         },
         metrics.adios.DeclareIO("EMFields")
      )
   );

   metrics.addMetric(
      std::make_unique<ParticleDumpMetric>(
         &g1,
         metrics.adios.DeclareIO(g1.name + "_dump")
      )
   );

   metrics.addMetric(
      std::make_unique<ParticleMetric>(
         &g1,
         metrics.adios.DeclareIO(g1.name + "_metrics"),
         Ncx, Ncy, Ncz
      )
   );

   return metrics;
}

template<int D>
double test_deposition(ParticleGroup& g1, const compute_t j_target) {
   std::println("Target J: {}", j_target);

   EMSolver emsolver(Nx, Ny, Nz, cfl, dt);
   constexpr BorisPush particle_push{};
   constexpr CurrentDeposition current_dep{};

   BorisPush::backstep_velocity(g1, emsolver.emdata);

   std::vector<compute_t> J(Nt);

   const auto metrics = create_metrics(
      "/home/cepheid/TriForce/game_engine/data/sp_jdep_test",
      emsolver, g1
   );

   metrics.write(0);
   for (std::size_t n = 0; n < Nt; n++) {
      emsolver.advance(dt * static_cast<compute_t>(n));
      particle_push(g1, emsolver.emdata, n);
      current_dep(emsolver.emdata, g1);

      metrics.write(n);

      g1.reset_y_positions();

      const auto& p = g1.particles[0];
      // std::println("{} {}", p.velocity, p.location);

      // if (n == 1) { exit(0); }

      assert(p.velocity[1] == 0.0_fp);
      assert(p.velocity[D] == 0.0_fp);
      const auto jx_sum = std::ranges::fold_left(emsolver.emdata.Jx.vec_data(), 0.0_fp, std::plus<compute_t>{});
      const auto jy_sum = std::ranges::fold_left(emsolver.emdata.Jy.vec_data(), 0.0_fp, std::plus<compute_t>{});
      const auto jz_sum = std::ranges::fold_left(emsolver.emdata.Jz.vec_data(), 0.0_fp, std::plus<compute_t>{});
      J[n] = jx_sum + jy_sum + jz_sum;
   }

   auto error = 0.0_fp;
   for (std::size_t i = 0; i < Nt; i++) { error += std::abs(J[i] - j_target); std::println("{}", J[i]); }
   return error / static_cast<compute_t>(Nt);
}

void electron_test_x() {
   const auto gamma = std::sqrt(1.0_fp + math::SQR(v_init / constants::c<compute_t>));

   std::println("Electron X-Deposition Test:");
   ParticleGroup g1("electrons", m_e, -q_e, 0);
   g1.particles.emplace_back(
      vec3{10.4875_fp, 0.5_fp, 25.5},  // location
      vec3{10.4875_fp, 0.5_fp, 25.5},  // old location
      vec3{v_init, 0.0_fp, 0.0_fp}, // velocity
      1.0_fp,                       // weight
      gamma,                        // gamma
      false                         // disabled
   );
   g1.initial_y_position = 0.5_fp;
   const auto j_target = J_target / gamma;
   const auto error = test_deposition<2>(g1, j_target);
   std::println("MAE = {:.3e} (≤ 4.340e-14)", error);
   // assert(error <= 4.340e-14);
}

void electron_test_z() {
   const auto gamma = std::sqrt(1.0_fp + math::SQR(v_init / constants::c<compute_t>));

   std::println("Electron Z-Deposition Test:");
   ParticleGroup g1("electrons", m_e, -q_e, 0);
   g1.particles.emplace_back(
      vec3{25.5_fp, 0.5_fp, 10.5},    // location
      vec3{25.5_fp, 0.5_fp, 10.5},  // old location
      vec3{0.0_fp, 0.0_fp, v_init}, // velocity
      1.0_fp,                         // weight
      gamma,                         // gamma
      false                           // disabled
   );
   g1.initial_y_position = 0.5_fp;
   const auto j_target = J_target / gamma;
   const auto error = test_deposition<0>(g1, j_target);
   std::println("MAE = {:.3e} (≤ 5.071e-14)", error);
   // assert(error <= 5.071e-14);
}

int main() {
   std::println("Single Particle Current Deposition Tests.");
   electron_test_x();
   // electron_test_z();
   return 0;
}