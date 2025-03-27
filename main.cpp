#include "program_params.hpp"
#include "constants.hpp"
#include "em_solver.hpp"
#include "metrics.hpp"
#include "sources.hpp"
#include "array.hpp"
#include "timers.hpp"
#include "particles.hpp"
#include "pusher.hpp"
#include "constants.hpp"
#include "current_deposition.hpp"

#include <print>
#include <chrono>

using tf::particles::ParticleGroup;
using tf::vec3;
using array_t = tf::Array3D<compute_t>;

// #include "particles.hpp"
// #include "octree.hpp"
// #include <vector>

int main() {
  //
  // std::vector<tf::particles::ParticleCell> cells((Nx - 1) * (Ny - 1) * (Nz - 1));
  //
  // for (std::size_t i = 0; i < Nx - 1; i++) {
  //   for (std::size_t j = 0; j < Ny - 1; j++) {
  //     for (std::size_t k = 0; k < Nz - 1; k++) {
  //       const auto code = tf::particles::morton_encode(i, j, k);
  //       cells[code] = {{}, code};
  //     }
  //   }
  // }
  // auto tree = tf::particles::create_particle_octree(cells);
  //
  // tf::particles::visit_octree(tree);

  // using continuous_t = tf::electromagnetics::ContinuousSource;
  // using temporal_vec = std::vector<std::unique_ptr<tf::electromagnetics::TemporalSource>>;

  // constexpr auto nx2 = Nx / 2;
  // constexpr auto freq = tf::constants::c / (20.0f * dx);


  // auto make_continuous = [&freq](temporal_vec& srcs) {
  //   srcs.push_back(std::make_unique<continuous_t>(2.0 * tf::constants::pi * freq, 0.0f, 0.0f, 1.0e30f, dx));
  // };

  // auto make_srcvec = [&]() -> temporal_vec {
  //   temporal_vec result{};
  //   // make_ricker(result);
  //   make_continuous(result);
  //   return result;
  // };

  const tf::electromagnetics::EMSolver emsolver(Nx, Ny, Nz, cfl, dt);

  // emsolver.emdata.srcs.emplace_back(
  //   &emsolver.emdata.Ez,
  //   tf::electromagnetics::SpatialSource{
  //     make_srcvec(),
  //     100.0f,
  //     {nx2, nx2 + 1, nx2, nx2 + 1, nx2, nx2 + 1}
  //   }
  // );
  //
  tf::metrics::Metrics metrics("/home/cepheid/TriForce/game_engine/data/particle_test");

  // metrics.addMetric(
  //   std::make_unique<tf::metrics::EMFieldsMetric>(
  //     std::unordered_map<std::string, array_t*> {
  //       {"Ex", &emsolver.emdata.Ex},
  //       {"Ey", &emsolver.emdata.Ey},
  //       {"Ez", &emsolver.emdata.Ez},
  //       {"Hx", &emsolver.emdata.Hx},
  //       {"Hy", &emsolver.emdata.Hy},
  //       {"Hz", &emsolver.emdata.Hz},
  //     },
  //     metrics.adios.DeclareIO("EMFields")
  //   )
  // );


  ParticleGroup g1("test_electron", tf::constants::m_e, -tf::constants::q_e, 0);

  tf::particles::ParticleInitializer::initializeFromFile(g1, "/home/cepheid/TriForce/game_engine/data/electron_init.txt");

  metrics.addMetric(
    std::make_unique<tf::metrics::ParticleMetric>(
      &g1,
      metrics.adios.DeclareIO("Particles")
    )
  );

  constexpr tf::particles::BorisPush particle_push{};

  Timer timer{};

  compute_t t = 0.0f;
  std::size_t step = 0zu;

  timer.start_timer();
  while (t <= total_time) {
    // emsolver.advance(t);

    particle_push(g1, emsolver.emdata);

    if (step % save_interval == 0) {
      const auto percent = 100.0f * t / total_time;
      std::println("Step {:4} Time: {:8.2e} Complete: {:4.1f}%", step, t, percent);
      metrics.write(step);
    }

    t += dt;
    step++;
  }
  const auto percent = 100.0f * t / total_time;
  std::println("Step {:4} Time: {:7.1e} Complete: {:4.1f}%", step, t, percent);
  metrics.write(step);

  timer.stop_timer();

  const auto time = std::chrono::hh_mm_ss(timer.elapsed);
  std::println("Total Time: {}", time);
  return 0;
}
