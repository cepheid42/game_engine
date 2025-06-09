#include "program_params.hpp"
#include "constants.hpp"
#include "em_solver.hpp"
#include "particles.hpp"
#include "pusher.hpp"
#include "vec3.hpp"
// #include "metrics.hpp"

#include <cassert>
#include <cmath>
#include <vector>
#include <algorithm>
#include <print>

using namespace tf;
using namespace tf::electromagnetics;
using namespace tf::particles;

// using namespace tf::metrics;
// Metrics create_metrics(const std::string& dir, const ParticleGroup& g1) {
//    Metrics metrics(dir);
//    metrics.addMetric(
//       std::make_unique<ParticleDumpMetric>(
//          &g1,
//          metrics.adios.DeclareIO(g1.name + "_dump")
//       )
//    );
//    return metrics;
// }


void electron_test() {
   constexpr auto m_e = constants::m_e<compute_t>;
   constexpr auto q_e = constants::q_e<compute_t>;
   ParticleGroup g1("electrons", m_e, -q_e, 0);

   constexpr auto init_loc = 11.5146959742731695_fp;
   constexpr auto init_vel =  1.186139936568e6_fp;

   g1.particles.emplace_back(
      vec3{7.5_fp, 0.5_fp, init_loc}, // location
      vec3{7.5_fp, 0.5_fp, init_loc}, // old location
      vec3{init_vel, 0.0_fp, 0.0_fp}, // velocity
      1.0_fp,                       // weight
      1.0_fp,                       // gamma
      false                         // disabled
   );

   g1.initial_y_position = 0.5_fp;

   EMSolver emsolver(Nx, Ny, Nz, cfl, dt);
   constexpr auto B_field = 4.2e-4_fp;
   emsolver.emdata.By_app.fill(B_field);
   emsolver.particle_correction();

   constexpr BorisPush particle_push{};
   BorisPush::backstep_velocity(g1, emsolver.emdata);

   std::vector<compute_t> xs(Nt);
   std::vector<compute_t> zs(Nt);

   for (std::size_t n = 0; n < Nt; n++) {
      particle_push(g1, emsolver.emdata, n);
      g1.reset_y_positions();
      xs[n] = x_range[0] + dx * g1.particles[0].location[0];
      zs[n] = z_range[0] + dz * g1.particles[0].location[2];
   }

   // for (int i = 0; i < Nt; i++) {
   //    std::println("[{}, {}],", xs[i], zs[i]);
   // }

   constexpr auto Te = 46424.0;
   const auto v_eth = std::sqrt(2.0 * constants::kB<compute_t> * Te / m_e);
   const auto gyro_radius = m_e * v_eth / (q_e * B_field);

   std::vector<compute_t> radius(Nt);
   std::ranges::transform(xs.begin(), xs.end(), zs.begin(), zs.end(),  radius.begin(),
      [&](const auto& x, const auto& z) { return std::sqrt(x * x + z * z); }
   );

   auto error = 0.0_fp;
   for (std::size_t i = 0; i < Nt; i++) {
      error += std::abs(gyro_radius - radius[i]);
   }
   error /= static_cast<compute_t>(Nt);
   std::println("Electron Gyroradius Test: MAE = {:.4e} <= 9.5e-9", error);
   assert(error <= 9.0e-9);
}

// void ion_test() {
//    constexpr auto m_p = constants::m_p<compute_t>;
//    constexpr auto q_p = constants::q_e<compute_t>;
//    ParticleGroup g1("ions", m_p, q_p, 1);
//    g1.particles.emplace_back(
//       vec3{2.5_fp, 0.5_fp, 5.5_fp}, // location
//       vec3{2.5_fp, 0.5_fp, 5.5_fp}, // old location
//       vec3{0.0_fp, 0.0_fp, 0.0_fp}, // velocity
//       1.0_fp,                       // weight
//       1.0_fp,                       // gamma
//       false                         // disabled
//    );
//
//    g1.initial_y_position = 0.5_fp;
//
//    EMSolver emsolver(Nx, Ny, Nz, cfl, dt);
//    constexpr auto E_field = 1.0_fp;
//    emsolver.emdata.Ex_app.fill(E_field);
//    emsolver.particle_correction();
//
//    constexpr BorisPush particle_push{};
//    BorisPush::backstep_velocity(g1, emsolver.emdata);
//
//    std::vector<compute_t> pos(Nt);
//
//    for (std::size_t n = 0; n < Nt; n++) {
//       particle_push(g1, emsolver.emdata, 1zu);
//       g1.reset_y_positions();
//
//       pos[n] = x_range[0] + dx * g1.particles[0].location[0];
//    }
//
//    const auto ts = math::linspace(0.0, t_end, Nt);
//    std::vector<compute_t> xs(Nt);
//    std::ranges::transform(ts, xs.begin(), [&](const auto& t) { return 0.5_fp * (q_p * E_field / m_p) * t * t; });
//
//    auto error = 0.0_fp;
//    for (int i = 0; i < Nt; i++) {
//       error += math::SQR(xs[i] - pos[i]);
//    }
//    error /= static_cast<compute_t>(Nt);
//    error = std::sqrt(error);
//
//    std::println("Ion Acceleration Test: RMSE = {:.3e} <= 2.0e-14", error);
//    assert(error <= 2.0e-14);
// }

int main() {
   std::println("Single Particle Gyro-Orbit Tests.");
   electron_test();
   // ion_test();
   return 0;
}