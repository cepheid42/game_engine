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
   BorisPush::backstep_velocity(g1, emsolver.emdata);
   std::vector<compute_t> pos(Nt);
   pos[0] = min_range[D] + deltas[D] * g1.particles[0].location[D];
   for (std::size_t n = 1; n < Nt; n++) {
      // no need to advance fields, since they don't change in time
      BorisPush::advance(g1, emsolver.emdata, 1zu);
      pos[n] = min_range[D] + deltas[D] * g1.particles[0].location[D];
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
   std::println("Ex Acceleration Test:\nMAE = {:.3e} (≤ 5.864e-11)", error);
   assert(error <= 5.864e-11);
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
   std::println("Ez Acceleration Test:\nMAE = {:.3e} (≤ 5.864e-11)", error);
   assert(error <= 5.864e-11);
}

int main() {
   std::println("Single Particle Acceleration Tests.");
   test_ex(-1.0, [](auto& emdata) { emdata.Ex_app.fill(-1.0); });
   test_ez(-1.0, [](auto& emdata) { emdata.Ez_app.fill(-1.0); });
   return 0;
}