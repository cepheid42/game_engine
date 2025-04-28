#include "program_params.hpp"
#include "constants.hpp"
#include "em_solver.hpp"
#include "metrics.hpp"
#include "array.hpp"
#include "timers.hpp"
#include "particles.hpp"
#include "pusher.hpp"
#include "current_deposition.hpp"

#include "barkeep.h"

#include <print>
#include <chrono>
#include <unordered_map>

using namespace tf;
using namespace tf::electromagnetics;
using namespace tf::particles;
using namespace tf::metrics;

namespace bk = barkeep;

void add_gaussianbeam(EMSolver& em) {
  using temporal_vec = std::vector<std::unique_ptr<TemporalSource>>;

  constexpr auto freq = static_cast<compute_t>(constants::c) / 8.0e-7_fp; // Hz -> c / 800 nm
  constexpr auto omega = 2.0_fp * static_cast<compute_t>(constants::pi) * freq;
  constexpr auto amp = 2.75e13_fp; // V/m
  constexpr auto w0 = 2.548e-6_fp; // meters, waste size

  constexpr auto width = 2.548e-14_fp; // seconds, ~25.48 fs
  constexpr auto delay = 2.0 * width;

  vec3 waist_pos{0.0_fp, 0.0_fp, 0.0_fp};

  constexpr auto x0 = 50zu;
  constexpr auto x1 = x0 + 1;
  constexpr auto y0 = 0zu;
  constexpr auto y1 = 1zu;
  constexpr auto z0 = 100zu;
  constexpr auto z1 = Nz - z0;

  using continuous_t = ContinuousSource;
  auto make_continuous = [&](temporal_vec& srcs) {
    srcs.push_back(std::make_unique<continuous_t>(omega, 0.0f, 0.0f, 1.0e30f, dx));
  };

  using gaussian_t = GaussianSource;
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
    SpatialSource(
      make_srcvec(),
      amp,
      {x0, x1, y0, y1, z0, z1}
    )
  );
}

Metrics create_metrics(const std::string& dir, EMSolver& em, const ParticleGroup& g1, const ParticleGroup& g2) {
  Metrics metrics(dir);
  
  metrics.addMetric(
    std::make_unique<EMFieldsMetric>(
      std::unordered_map<std::string, Array3D<compute_t>*> {
        {"Ex", &em.emdata.Ex},
        {"Ey", &em.emdata.Ey},
        {"Ez", &em.emdata.Ez},
        {"Hx", &em.emdata.Hx},
        {"Hy", &em.emdata.Hy},
        {"Hz", &em.emdata.Hz},
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
      metrics.adios.DeclareIO(g1.name + "_metrics")
    )
  );

  metrics.addMetric(
    std::make_unique<ParticleDumpMetric>(
      &g2,
      metrics.adios.DeclareIO(g2.name + "_dump")
    )
  );

  metrics.addMetric(
    std::make_unique<ParticleMetric>(
      &g2,
      metrics.adios.DeclareIO(g2.name + "_metrics")
    )
  );

  return metrics;
}

std::unordered_map<std::string, Timer> create_timers() {
  std::unordered_map<std::string, Timer> timers{};
  timers["Main"] = Timer{};
  timers["EM"] = Timer{};
  timers["Push"] = Timer{};
  timers["Jdep"] = Timer{};
  timers["IO"] = Timer{};
  return timers;
}

void print_final_timers(auto& timers) {
  std::println("   EM: {}", std::chrono::hh_mm_ss(timers["EM"].elapsed));
  std::println(" Push: {}", std::chrono::hh_mm_ss(timers["Push"].elapsed));
  std::println(" Jdep: {}", std::chrono::hh_mm_ss(timers["Jdep"].elapsed));
  std::println("   IO: {}", std::chrono::hh_mm_ss(timers["IO"].elapsed));
  std::println("Total: {}", std::chrono::hh_mm_ss(timers["Main"].elapsed));
}


int main() {
  auto timers = create_timers();

  timers["Main"].start_timer();
  constexpr auto m_e = static_cast<compute_t>(constants::m_e);
  constexpr auto m_p = static_cast<compute_t>(constants::m_p);
  constexpr auto q_e = static_cast<compute_t>(constants::q_e);
  constexpr auto ion_file = "/home/cepheid/TriForce/game_engine/data/ion_slab.dat";
  constexpr auto electron_file = "/home/cepheid/TriForce/game_engine/data/electron_slab.dat";

  auto g1 = ParticleInitializer::initializeFromFile("electrons", m_e, -q_e, 0, electron_file);
  auto g2 = ParticleInitializer::initializeFromFile("ions", m_p, +q_e, 1, ion_file);

  EMSolver emsolver(Nx, Ny, Nz, cfl, dt);
  add_gaussianbeam(emsolver);
  const CurrentDeposition jdep{emsolver.emdata};
  constexpr BorisPush particle_push{};

  auto metrics = create_metrics("/home/cepheid/TriForce/game_engine/data/lsi_test", emsolver, g1, g2);

  compute_t t = 0.0_fp;
  std::size_t step = 0zu;

  const auto progress_bar =
    bk::ProgressBar(&step, {
      .total = Nt,
      .message = "Step",
      .speed = 0.,
      .speed_unit = "steps/s",
      .interval = 1.,
      // .no_tty = true,
      .show = false}
    );

  timers["IO"].start_timer();
  metrics.write(step);
  timers["IO"].stop_timer();

  progress_bar->show();
  while (t <= total_time) {
    timers["EM"].start_timer();
    emsolver.advance(t);
    timers["EM"].stop_timer();

    timers["Push"].start_timer();
    particle_push(g1, emsolver.emdata);
    particle_push(g2, emsolver.emdata);
    timers["Push"].stop_timer();

    timers["Jdep"].start_timer();
    jdep(g1, emsolver.emdata);
    jdep(g2, emsolver.emdata);
    timers["Jdep"].stop_timer();

    g1.reset_y_positions();
    g2.reset_y_positions();

    t += dt;
    step++;

    if (step % save_interval == 0) {
      timers["IO"].start_timer();
      metrics.write(step);
      timers["IO"].stop_timer();
    }
  }
  progress_bar->done();
  timers["Main"].stop_timer();

  print_final_timers(timers);
  return 0;
}
