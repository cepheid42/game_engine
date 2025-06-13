#include "compute_type.hpp"
#include "program_params.hpp"
#include "constants.hpp"
#include "electromagnetics/em_solver.hpp"
#include "metrics.hpp"
#include "array.hpp"
#include "timers.hpp"
#include "particles/particles.hpp"
#include "particles/pusher.hpp"
#include "particles/current_deposition.hpp"

#include "barkeep.h"

// #include <print>

using namespace tf;
using namespace tf::electromagnetics;
using namespace tf::particles;
using namespace tf::metrics;

namespace bk = barkeep;

Metrics create_metrics(const std::string& dir, EMSolver& em, const ParticleGroup& g1, const ParticleGroup& g2) {
   Metrics metrics(dir);

   metrics.addMetric(
      std::make_unique<EMFieldsMetric>(
         std::unordered_map<std::string, Array3D<compute_t>*>{
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
         metrics.adios.DeclareIO(g1.name + "_metrics"),
         Ncx, Ncy, Ncz
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
         metrics.adios.DeclareIO(g2.name + "_metrics"),
         Ncx, Ncy, Ncz
      )
   );

   return metrics;
}

int main() {
   auto timers = utilities::create_timers();

   timers["Main"].start_timer();
   constexpr auto m_e = constants::m_e<compute_t>;
   constexpr auto m_p = constants::m_p<compute_t>;
   constexpr auto q_e = constants::q_e<compute_t>;
   constexpr auto electron_file = "/home/cepheid/TriForce/game_engine/data/electron_slab.dat";
   constexpr auto ion_file = "/home/cepheid/TriForce/game_engine/data/ion_slab.dat";

   auto g1 = ParticleInitializer::initializeFromFile("electrons", m_e, -q_e, 0, electron_file);
   auto g2 = ParticleInitializer::initializeFromFile("ions", m_p, +q_e, 1, ion_file);

   EMSolver emsolver(Nx, Ny, Nz, cfl, dt);
   add_gaussianbeam(emsolver);

   constexpr BorisPush particle_push{};
   constexpr CurrentDeposition current_dep{};
   BorisPush::backstep_velocity(g1, emsolver.emdata);
   BorisPush::backstep_velocity(g2, emsolver.emdata);

   const auto metrics = create_metrics(
      "/home/cepheid/TriForce/game_engine/data/lsi_test",
      emsolver,
      g1, g2
   );

   std::size_t step = 0zu;
   const auto progress_bar =
      bk::ProgressBar(
         &step,
         {.total = Nt,
          .message = "Step",
          .speed = 0.,
          .speed_unit = "steps/s",
          .interval = 1.,
          .no_tty = true,
          .show = false});

   timers["IO"].start_timer();
   metrics.write(step);
   timers["IO"].stop_timer();

   progress_bar->show();
   for (step = 0; step < Nt; step++) {
      timers["EM"].start_timer();
      emsolver.advance(static_cast<compute_t>(step) * dt);
      timers["EM"].stop_timer();

      timers["Push"].start_timer();
      particle_push(g1, emsolver.emdata, step);
      particle_push(g2, emsolver.emdata, step);
      timers["Push"].stop_timer();

      timers["Jdep"].start_timer();
      current_dep(emsolver.emdata, g1);
      current_dep(emsolver.emdata, g2);
      timers["Jdep"].stop_timer();

      g1.reset_y_positions();
      g2.reset_y_positions();

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
