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
#include <fstream>

using namespace tf;
using namespace tf::electromagnetics;
using namespace tf::particles;

inline constexpr auto m_e = constants::m_e<compute_t>;
inline constexpr auto q_e = constants::q_e<compute_t>;

void electron_test() {
   ParticleGroup g1("electrons", m_e, -q_e, 0);

   constexpr auto init_loc = 11.5146959742731695_fp;
   constexpr auto init_vel =  1.186139936568e6_fp;
   constexpr auto B_field = 4.2e-4_fp;

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
   emsolver.emdata.By_app.fill(B_field);
   emsolver.particle_correction();
   BorisPush::backstep_velocity(g1, emsolver.emdata);
   std::vector<compute_t> xs(Nt);
   std::vector<compute_t> zs(Nt);

   for (std::size_t n = 0; n < Nt; n++) {
      BorisPush::advance(g1, emsolver.emdata, n);
      xs[n] = x_range[0] + dx * g1.particles[0].location[0];
      zs[n] = z_range[0] + dz * g1.particles[0].location[2];
   }

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
   std::println("Electron Gyroradius Test:\nMAE = {:.3e} (≤ 8.246e-9)", error);
   assert(error <= 8.246e-9);
}

void larmor_test() {
   // const auto xr = math::linspace(x_range[0], x_range[1], Nx);
   // const auto zr = math::linspace(z_range[0], z_range[1], Nz);
   //
   // std::vector<compute_t> Ex{};
   // std::vector<compute_t> Ez{};
   // for (int i = 0; i < Nx - 1; i++) {
   //    for (int k = 0; k < Nz; k++) {
   //       const auto r = std::sqrt(math::SQR(xr[i]) + math::SQR(zr[k] + 0.5 * dz));
   //       const auto theta = std::atan2(zr[k] + 0.5 * dz, xr[i]);
   //       Ex.push_back((15 / r) * std::cos(theta));
   //    }
   // }
   //
   // for (int i = 0; i < Nx; i++) {
   //    for (int k = 0; k < Nz - 1; k++) {
   //       const auto r = std::sqrt(math::SQR(xr[i] + 0.5 * dx) + math::SQR(zr[k]));
   //       const auto theta = std::atan2(zr[k], xr[i] + 0.5 * dz);
   //       Ez.push_back((15 / r) * std::sin(theta));
   //    }
   // }

   ParticleGroup g1("electrons", m_e, -q_e, 0);
   constexpr auto B_field = 7.0e-4_fp;
   // constexpr auto init_loc = 4.5_fp;
   constexpr auto init_vel = -3.0e5_fp;
   g1.particles.emplace_back(
      vec3{32.5_fp, 0.5_fp, 50.5}, // location
      vec3{32.5_fp, 0.5_fp, 50.5}, // old location
      vec3{init_vel, 0.0_fp, 0.0_fp}, // velocity
      1.0_fp,                       // weight
      1.0_fp,                       // gamma
      false                         // disabled
   );

   g1.initial_y_position = 0.5_fp;

   EMSolver emsolver(Nx, Ny, Nz, cfl, dt);
   // std::ranges::copy(Ex, emsolver.emdata.Ex_app.data());
   // std::ranges::copy(Ez, emsolver.emdata.Ez_app.data());
   emsolver.emdata.By_app.fill(B_field);
   emsolver.emdata.Ex_app.fill(-25);
   emsolver.emdata.Ez_app.fill(-25);
   emsolver.particle_correction();
   BorisPush::backstep_velocity(g1, emsolver.emdata);

   std::vector<compute_t> xs(Nt);
   std::vector<compute_t> zs(Nt);

   for (std::size_t n = 0; n < Nt; n++) {
      BorisPush::advance(g1, emsolver.emdata, n);
      xs[n] = x_range[0] + dx * g1.particles[0].location[0];
      zs[n] = z_range[0] + dz * g1.particles[0].location[2];
   }

   // constexpr auto Te = 46424.0;
   // const auto v_eth = std::sqrt(2.0 * constants::kB<compute_t> * Te / m_e);
   // const auto gyro_radius = m_e * v_eth / (q_e * B_field);
   //
   std::vector<compute_t> radius(Nt);
   std::ranges::transform(xs.begin(), xs.end(), zs.begin(), zs.end(),  radius.begin(),
      [&](const auto& x, const auto& z) { return std::sqrt(x * x + z * z); }
   );

   std::println("Inner = {}, Outer = {}", std::ranges::min(radius), std::ranges::max(radius));

   std::ofstream file("/home/cepheid/TriForce/game_engine/data/larmor.csv");
   file << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10);
   for (int i = 0; i < xs.size(); i++) {
      file << xs[i] << ", " << zs[i] << std::endl;
   }
   // for (int i = 0; i < Ex.size(); i++) {
   //    file << Ex[i] << ", " << Ez[i] << std::endl;
   // }
   file.close();

   // for (const auto& x : xs) {
   //    std::print("{}, ", x);
   // }
   // std::println();
   // for (const auto& x : zs) {
   //    std::print("{}, ", x);
   // }

   // auto error = 0.0_fp;
   // for (std::size_t i = 0; i < Nt; i++) {
   //    error += std::abs(gyro_radius - radius[i]);
   // }
   // error /= static_cast<compute_t>(Nt);
   // std::println("Electron Gyroradius Test:\nMAE = {:.3e} (≤ 8.246e-9)", error);
   // assert(error <= 8.246e-9);
}


int main() {
   std::println("Single Particle Gyro-Orbit Tests.");
   // electron_test();
   larmor_test();
   return 0;
}