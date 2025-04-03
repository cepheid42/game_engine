#include "program_params.hpp"
#include "constants.hpp"
#include "em_solver.hpp"
#include "metrics.hpp"
#include "array.hpp"
#include "timers.hpp"
#include "particles.hpp"
#include "pusher.hpp"
// #include "current_deposition.hpp"
// #include "constants.hpp"

#include <print>
#include <chrono>

using tf::particles::ParticleGroup;
using tf::vec3;
using array_t = tf::Array3D<compute_t>;

void add_em_sources(tf::electromagnetics::EMSolver& em) {
  using temporal_vec = std::vector<std::unique_ptr<tf::electromagnetics::TemporalSource>>;

  constexpr auto nx2 = Nx / 2;
  constexpr auto freq = tf::constants::c / (20.0f * dx);

  using ricker_t = tf::electromagnetics::RickerSource;
  auto make_ricker = [&freq](temporal_vec& srcs) {
    srcs.push_back(std::make_unique<ricker_t>(freq));
  };

  // using continuous_t = tf::electromagnetics::ContinuousSource;
  // auto make_continuous = [&freq](temporal_vec& srcs) {
  //   srcs.push_back(std::make_unique<continuous_t>(2.0 * tf::constants::pi * freq, 0.0f, 0.0f, 1.0e30f, dx));
  // };

  auto make_srcvec = [&]() -> temporal_vec {
    temporal_vec result{};
    make_ricker(result);
    // make_continuous(result);
    return result;
  };

  em.emdata.srcs.emplace_back(
    &em.emdata.Ez,
    tf::electromagnetics::SpatialSource{
      make_srcvec(),
      100.0f,
      {nx2 + 14, nx2 + 15, nx2 + 14, nx2 + 15, nx2, nx2 + 1}
    }
  );
}


int main() {

  tf::electromagnetics::EMSolver emsolver(Nx, Ny, Nz, cfl, dt);

  // add_em_sources(emsolver);

  // for (auto& e: emsolver.emdata.Bz_app) {
  //   e = 1600.0e-4f;
  // }

  // for (auto& e: emsolver.emdata.Ez) {
  //   e = 1000.0f;
  // }

  tf::metrics::Metrics metrics("/home/cepheid/TriForce/game_engine/data/particle_test");

  // metrics.addMetric(
  //   std::make_unique<tf::metrics::EMFieldsMetric>(
  //     std::unordered_map<std::string, array_t*> {
  //       // {"Ex", &emsolver.emdata.Ex},
  //       // {"Ey", &emsolver.emdata.Ey},
  //       {"Ez", &emsolver.emdata.Ez},
  //       // {"Hx", &emsolver.emdata.Hx},
  //       // {"Hy", &emsolver.emdata.Hy},
  //       // {"Hz", &emsolver.emdata.Hz},
  //     },
  //     metrics.adios.DeclareIO("EMFields")
  //   )
  // );

  using tf::particles::ParticleInitializer;
  constexpr auto m_e = static_cast<compute_t>(tf::constants::m_e);
  constexpr auto q_e = static_cast<compute_t>(tf::constants::q_e);
  auto g1 = ParticleInitializer::initializeFromFile("electrons", m_e, -q_e, 0, "/home/cepheid/TriForce/game_engine/data/electrons.dat");

  metrics.addMetric(
    std::make_unique<tf::metrics::ParticleMetric>(
      &g1,
      metrics.adios.DeclareIO("Particles")
    )
  );

  constexpr tf::particles::BorisPush particle_push{};
  // constexpr tf::particles::CurrentDeposition jdep{};

  Timer timer{};

  compute_t t = 0.0_fp;
  std::size_t step = 0zu;

  timer.start_timer();

  while (t <= total_time) {
    emsolver.advance(t);

    particle_push(g1, emsolver.emdata);

    // jdep(g1, emsolver.emdata);

    if (step % save_interval == 0) {
      const auto percent = 100.0_fp * t / total_time;
      std::println("Step {:4} Time: {:8.2e} Complete: {:4.1f}%", step, t, percent);
      metrics.write(step);
    }

    t += dt;
    step++;
  }
  const auto percent = 100.0_fp * t / total_time;
  std::println("Step {:4} Time: {:7.1e} Complete: {:4.1f}%", step, t, percent);
  metrics.write(step);

  timer.stop_timer();

  const auto time = std::chrono::hh_mm_ss(timer.elapsed);
  std::println("Total Time: {}", time);
  return 0;
}
