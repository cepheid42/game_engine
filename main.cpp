#include "program_params.hpp"
#include "constants.hpp"
#include "em_solver.hpp"
#include "metrics.hpp"
#include "array.hpp"
#include "timers.hpp"
#include "particles.hpp"
#include "pusher.hpp"
#include "current_deposition.hpp"

// #define UTL_PROFILER_DISABLE
#include "profiler.hpp"

#include <print>
#include <chrono>

using tf::particles::ParticleGroup;
using tf::vec3;
using array_t = tf::Array3D<compute_t>;

void add_gaussianbeam(tf::electromagnetics::EMSolver& em) {
  using temporal_vec = std::vector<std::unique_ptr<tf::electromagnetics::TemporalSource>>;

  constexpr auto freq = static_cast<compute_t>(tf::constants::c) / 8.0e-7_fp; // Hz -> c / 800 nm
  constexpr auto omega = 2.0_fp * static_cast<compute_t>(tf::constants::pi) * freq;
  constexpr auto amp = 2.75e13_fp; // V/m
  constexpr auto w0 = 2.548e-6_fp; // meters, waste size

  constexpr auto width = 2.548e-14_fp; // seconds, ~25.48 fs
  constexpr auto delay = 2.0 * width;

  vec3 waist_pos{0.0_fp, 0.0_fp, 0.0_fp};

  constexpr auto x0 = 100zu;
  constexpr auto x1 = x0 + 1;
  constexpr auto y0 = 0zu;
  constexpr auto y1 = 1zu;
  constexpr auto z0 = 100zu;
  constexpr auto z1 = Nz - z0;

  using continuous_t = tf::electromagnetics::ContinuousSource;
  auto make_continuous = [&](temporal_vec& srcs) {
    srcs.push_back(std::make_unique<continuous_t>(omega, 0.0f, 0.0f, 1.0e30f, dx));
  };

  using gaussian_t = tf::electromagnetics::GaussianSource;
  auto make_gaussian = [&](temporal_vec& srcs) {
    srcs.push_back(std::make_unique<gaussian_t>(width, 2.0_fp, delay));
  };


  auto make_srcvec = [&]() -> temporal_vec {
    temporal_vec result{};
    make_gaussian(result);
    make_continuous(result);
    return result;
  };

  em.emdata.srcs.emplace_back(
    &em.emdata.Ey,
    w0,
    omega,
    waist_pos,
    tf::electromagnetics::SpatialSource(
      make_srcvec(),
      amp,
      {x0, x1, y0, y1, z0, z1}
    )
  );
}


int main() {
  UTL_PROFILER_SCOPE("Main");
  // std::println("dx, dy, dz: {}, {}, {}", dx, dy, dz);
  // std::println("dt: {}", dt);
  // std::println("Total time: {}", total_time);
  // std::println("Nt: {}", static_cast<int>(total_time / dt));

  tf::electromagnetics::EMSolver emsolver(Nx, Ny, Nz, cfl, dt);

  add_gaussianbeam(emsolver);

  tf::metrics::Metrics metrics("/home/cepheid/TriForce/game_engine/data/lsi_test");

  metrics.addMetric(
    std::make_unique<tf::metrics::EMFieldsMetric>(
      std::unordered_map<std::string, array_t*> {
        {"Ex", &emsolver.emdata.Ex},
        {"Ey", &emsolver.emdata.Ey},
        {"Ez", &emsolver.emdata.Ez},
        {"Hx", &emsolver.emdata.Hx},
        {"Hy", &emsolver.emdata.Hy},
        {"Hz", &emsolver.emdata.Hz},
      },
      metrics.adios.DeclareIO("EMFields")
    )
  );

  // using tf::particles::ParticleInitializer;
  // constexpr auto m_e = static_cast<compute_t>(tf::constants::m_e);
  // constexpr auto q_e = static_cast<compute_t>(tf::constants::q_e);
  // constexpr auto ion_file = "/home/cepheid/TriForce/game_engine/data/ion_slab.dat";
  // constexpr auto electron_file = "/home/cepheid/TriForce/game_engine/data/electron_slab.dat";
  //
  //
  // auto g1 = ParticleInitializer::initializeFromFile("electrons", m_e, -q_e, 0, electron_file);
  // auto g2 = ParticleInitializer::initializeFromFile("ions", 1836.152674_fp * m_e, q_e, 1, ion_file);
  //
  // metrics.addMetric(
  //   std::make_unique<tf::metrics::ParticleDumpMetric>(
  //     &g1,
  //     metrics.adios.DeclareIO(g1.name + "_dump")
  //   )
  // );
  //
  // metrics.addMetric(
  //   std::make_unique<tf::metrics::ParticleMetric>(
  //     &g1,
  //     metrics.adios.DeclareIO(g1.name + "_metrics")
  //   )
  // );
  //
  // metrics.addMetric(
  //   std::make_unique<tf::metrics::ParticleDumpMetric>(
  //     &g2,
  //     metrics.adios.DeclareIO(g2.name + "_dump")
  //   )
  // );
  //
  // metrics.addMetric(
  //   std::make_unique<tf::metrics::ParticleMetric>(
  //     &g2,
  //     metrics.adios.DeclareIO(g2.name + "_metrics")
  //   )
  // );
  //
  // constexpr tf::particles::BorisPush particle_push{};
  // tf::particles::CurrentDeposition jdep{emsolver.emdata};

  Timer main_timer{};
  Timer em_timer{};
  Timer push_timer{};
  Timer jdep_timer{};
  Timer metrics_timer{};

  compute_t t = 0.0_fp;
  std::size_t step = 0zu;

  main_timer.start_timer();

  std::println("Step {:4} Time: {:8.2e} Complete: {:4.1f}%", step, t, 0.0_fp);
  metrics_timer.start_timer();
  metrics.write(step);
  metrics_timer.stop_timer();

  // constexpr auto freq = static_cast<compute_t>(tf::constants::c) / (20.0_fp * dx);
  // const tf::electromagnetics::RickerSource ricker{freq};

  while (t <= total_time) {
    UTL_PROFILER_BEGIN(em, "Electromagnetics");
    em_timer.start_timer();
    emsolver.advance(t);
    em_timer.stop_timer();
    UTL_PROFILER_END(em);

    // UTL_PROFILER_BEGIN(push, "Particle Push");
    // push_timer.start_timer();
    // particle_push(g1, emsolver.emdata);
    // particle_push(g2, emsolver.emdata);
    // push_timer.stop_timer();
    // UTL_PROFILER_END(push);
    //
    // UTL_PROFILER_BEGIN(jdep, "Particle Deposition");
    // jdep_timer.start_timer();
    // jdep(g1, emsolver.emdata);
    // jdep(g2, emsolver.emdata);
    // jdep_timer.stop_timer();
    // UTL_PROFILER_END(jdep);

    t += dt;
    step++;

    if (step % save_interval == 0) {
      metrics_timer.start_timer();
      const auto percent = 100.0_fp * t / total_time;
      std::println("Step {:4} Time: {:8.2e} Complete: {:4.1f}%", step, t, percent);
      metrics.write(step);
      metrics_timer.stop_timer();
    }
  }

  main_timer.stop_timer();

  std::println("   EM: {}", std::chrono::hh_mm_ss(em_timer.elapsed));
  std::println(" Push: {}", std::chrono::hh_mm_ss(push_timer.elapsed));
  std::println(" Jdep: {}", std::chrono::hh_mm_ss(jdep_timer.elapsed));
  std::println("   IO: {}", std::chrono::hh_mm_ss(metrics_timer.elapsed));
  std::println("Total: {}", std::chrono::hh_mm_ss(main_timer.elapsed));
  return 0;
}
