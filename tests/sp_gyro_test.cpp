#include "program_params.hpp"
#include "constants.hpp"
#include "em_solver.hpp"
#include "particles.hpp"
#include "pusher.hpp"
#include "current_deposition.hpp"
#include "vec3.hpp"

#include <cassert>
#include <cmath>
#include <vector>
#include <algorithm>
#include <print>
#include <iomanip>

using namespace tf;
using namespace tf::electromagnetics;
using namespace tf::particles;

inline constexpr auto m_e = constants::m_e<compute_t>;
inline constexpr auto q_e = constants::q_e<compute_t>;

void location2File(const auto& x, const auto& z, const std::string& path) {
   std::ofstream file(path);
   file << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10);
   for (int i = 0; i < Nt; i++) {
      file << x[i] << ", " << 0.0 << ", " << z[i] << std::endl;
   }
   file.close();
}

auto gyroradius_test() {
   constexpr auto B_field = 4.2e-4_fp;
   constexpr auto Te = 46424.0;
   const auto v_eth = vec3{std::sqrt(2.0 * constants::kB<compute_t> * Te / m_e), 0.0, 0.0};
   const auto gamma = calculateGamma(v_eth);
   const auto gyro_radius = gamma * m_e * v_eth.length() / (q_e * B_field);
   constexpr auto init_loc = 11.51475147194894_fp;

   ParticleGroup g1("electrons", m_e, -q_e, 0);

   g1.particles.emplace_back(
      vec3{7.5_fp, 0.5_fp, init_loc}, // location
      vec3{7.5_fp, 0.5_fp, init_loc}, // old location
      v_eth,                          // velocity
      1.0_fp,                       // weight
      1.0_fp,                       // gamma
      false                         // disabled
   );

   g1.initial_y_position = 0.5_fp;

   EMSolver emsolver(Nx, Ny, Nz, cfl, dt);
   emsolver.emdata.By_app.fill(B_field);
   emsolver.particle_correction();
   BorisPush::backstep_velocity(g1, emsolver.emdata);
   std::vector<compute_t> xp(Nt);
   std::vector<compute_t> zp(Nt);

   for (std::size_t n = 0; n < Nt; n++) {
      xp[n] = x_range[0] + dx * g1.particles[0].location[0];
      zp[n] = z_range[0] + dz * g1.particles[0].location[2];
      BorisPush::advance(g1, emsolver.emdata, n);
   }

   std::vector<compute_t> radius(Nt);
   std::ranges::transform(xp.begin(), xp.end(), zp.begin(), zp.end(),  radius.begin(),
      [&](const auto& x, const auto& z) { return std::sqrt(x * x + z * z); }
   );

   auto error = 0.0_fp;
   for (std::size_t i = 0; i < Nt; i++) {
      error += std::abs(gyro_radius - radius[i]);
   }
   error /= static_cast<compute_t>(Nt);
   std::println("Electron Gyroradius Test:\n\tMAE = {:.3e} (≤ 7.039e-8)", error);
   assert(error <= 7.039e-8);
}

void ExBy_test() {
   constexpr auto B_field = 2.0e-3_fp;
   constexpr auto E_field = -150.0;
   constexpr auto v_init = 1.0e6_fp;

   ParticleGroup g1("electrons", m_e, -q_e, 0);
   g1.particles.emplace_back(
      vec3{7.5_fp, 0.5_fp, 13.5},     // location
      vec3{7.5_fp, 0.5_fp, 13.5},     // old location
      vec3{v_init, 0.0_fp, 0.0_fp}, // velocity
      1.0_fp,                         // weight
      1.0_fp,                         // gamma
      false                           // disabled
   );

   g1.initial_y_position = 0.5_fp;

   EMSolver emsolver(Nx, Ny, Nz, cfl, dt);
   emsolver.emdata.By_app.fill(B_field);
   emsolver.emdata.Ex_app.fill(E_field);
   emsolver.particle_correction();
   BorisPush::backstep_velocity(g1, emsolver.emdata);

   std::vector<compute_t> xp(Nt);
   std::vector<compute_t> zp(Nt);
   xp[0] = x_range[0] + dx * g1.particles[0].location[0];
   zp[0] = z_range[0] + dz * g1.particles[0].location[2];
   for (std::size_t n = 1; n < Nt; n++) {
      BorisPush::advance(g1, emsolver.emdata, n);
      xp[n] = x_range[0] + dx * g1.particles[0].location[0];
      zp[n] = z_range[0] + dz * g1.particles[0].location[2];
   }

   // location2File(xp, zp, "/home/cepheid/TriForce/game_engine/data/larmor.csv");

   std::println("Ex/By ExB Drift");
   constexpr auto gyro_radius = m_e * v_init / (q_e * B_field);
   const auto xmin = std::ranges::min(xp);
   const auto xmax = std::ranges::max(xp);
   const auto lr = (xmax - xmin) / 2.0;
   std::println("\tGyroradius error = {:.3e} (≤ 8.426.e-6)", std::abs(lr - gyro_radius));
   assert(std::abs(lr - gyro_radius) <= 8.426e-6);

   const auto zmin = std::ranges::min(zp);
   const auto zmax = std::ranges::max(zp);
   constexpr auto v_gc = E_field / B_field;
   const auto z_end = std::abs(v_gc * t_end);
   const auto deltaZ = (zmax - zmin) - 2.0 * gyro_radius;
   std::println("\tDrift error = {:.3e} (≤ 5.887-6)", std::abs(z_end - deltaZ));
   assert(std::abs(z_end - deltaZ) <= 5.887e-6);
}

void EzBy_test() {
   constexpr auto B_field = 2.0e-3_fp;
   constexpr auto E_field = -150.0;
   constexpr auto v_init = 1.0e6_fp;

   ParticleGroup g1("electrons", m_e, -q_e, 0);
   g1.particles.emplace_back(
      vec3{1.5_fp, 0.5_fp, 7.5},     // location
      vec3{1.5_fp, 0.5_fp, 7.5},     // old location
      vec3{0.0_fp, 0.0_fp, v_init}, // velocity
      1.0_fp,                         // weight
      1.0_fp,                         // gamma
      false                           // disabled
   );

   g1.initial_y_position = 0.5_fp;

   EMSolver emsolver(Nx, Ny, Nz, cfl, dt);
   emsolver.emdata.By_app.fill(B_field);
   emsolver.emdata.Ez_app.fill(E_field);
   emsolver.particle_correction();
   BorisPush::backstep_velocity(g1, emsolver.emdata);

   std::vector<compute_t> xp(Nt);
   std::vector<compute_t> zp(Nt);
   xp[0] = x_range[0] + dx * g1.particles[0].location[0];
   zp[0] = z_range[0] + dz * g1.particles[0].location[2];
   for (std::size_t n = 1; n < Nt; n++) {
      BorisPush::advance(g1, emsolver.emdata, n);

      xp[n] = x_range[0] + dx * g1.particles[0].location[0];
      zp[n] = z_range[0] + dz * g1.particles[0].location[2];
   }

   // location2File(xp, zp, "/home/cepheid/TriForce/game_engine/data/larmor.csv");

   std::println("Ez/By ExB Drift");
   constexpr auto gyro_radius = m_e * v_init / (q_e * B_field);
   const auto zmin = std::ranges::min(zp);
   const auto zmax = std::ranges::max(zp);
   const auto lr = (zmax - zmin) / 2.0;
   std::println("\tGyroradius error = {:.3e} (≤ 8.426.e-6)", std::abs(lr - gyro_radius));
   assert(std::abs(lr - gyro_radius) <= 8.426e-6);

   const auto xmin = std::ranges::min(xp);
   const auto xmax = std::ranges::max(xp);
   constexpr auto v_gc = E_field / B_field;
   const auto x_end = std::abs(v_gc * t_end);
   const auto deltaX = (xmax - xmin) - 2.0 * gyro_radius;
   std::println("\tDrift error = {:.3e} (≤ 5.887-6)", std::abs(x_end - deltaX));
   assert(std::abs(x_end - deltaX) <= 5.887e-6);
}

int main() {
   std::println("Single Particle Gyro-Orbit Tests.");
   gyroradius_test();
   ExBy_test();
   EzBy_test();
   return 0;
}