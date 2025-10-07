#include "program_params.hpp"
#include "em_solver.hpp"
#include "metrics.hpp"
#include "timers.hpp"
#include "particles/particles.hpp"
#include "particles/pusher.hpp"
#include "particles/current_deposition.hpp"

#include "barkeep.h"

using namespace tf;
using namespace tf::electromagnetics;
using namespace tf::particles;
using namespace tf::metrics;

namespace bk = barkeep;

void add_em_metrics(Metrics& metrics, auto& em) {
   metrics.addMetric(
      std::make_unique<EMFieldsMetric>(
         std::unordered_map<std::string, Array3D<double>*>{
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
}

void add_group_metric(Metrics& metrics, const auto& pg) {
   metrics.addMetric(
      std::make_unique<ParticleDumpMetric>(
         &pg,
         metrics.adios.DeclareIO(pg.name + "_dump")
      )
   );

   // metrics.addMetric(
   //    std::make_unique<ParticleMetric>(
   //       &pg,
   //       metrics.adios.DeclareIO(pg.name + "_metrics"),
   //       Nx - 1, Ny - 1, Nz - 1
   //    )
   // );
}


int main() {
   auto timers = utilities::create_timers();
   timers["Main"].start_timer();
   constexpr auto electron_file = "/home/cepheid/TriForce/game_engine/data/electrons.dat";
   constexpr auto      ion_file = "/home/cepheid/TriForce/game_engine/data/ions.dat";
   auto g1 = ParticleInitializer::initializeFromFile(electron_file);
   auto g2 = ParticleInitializer::initializeFromFile(ion_file);

   EMSolver emsolver(Nx, Ny, Nz);
   add_gaussianbeam(emsolver);

   emsolver.particle_correction();
   BorisPush::backstep_velocity(g1, emsolver.emdata);
   BorisPush::backstep_velocity(g2, emsolver.emdata);

   Metrics metrics("/home/cepheid/TriForce/game_engine/data/lsi_test");
   add_em_metrics(metrics, emsolver);
   add_group_metric(metrics, g1);
   add_group_metric(metrics, g2);

   auto step = 0lu;
   const auto progress_bar =
      bk::ProgressBar(
         &step,
         {.total = Nt,
          .message = "Step",
          .speed = 0.,
          .speed_unit = "steps/s",
          .interval = 1.,
          // .no_tty = true,
          .show = false});

   timers["IO"].start_timer();
   metrics.write(step);
   timers["IO"].stop_timer();

   progress_bar->show();
   for (step = 1; step <= Nt; step++) {
      timers["EM"].start_timer();
      emsolver.advance(static_cast<double>(step) * dt);
      timers["EM"].stop_timer();

      timers["Push"].start_timer();
      g1.reset_y_positions();
      g2.reset_y_positions();
      BorisPush::advance(g1, emsolver.emdata, step);
      BorisPush::advance(g2, emsolver.emdata, step);
      timers["Push"].stop_timer();

      timers["Jdep"].start_timer();
      CurrentDeposition::advance(g1, emsolver.emdata);
      CurrentDeposition::advance(g2, emsolver.emdata);
      timers["Jdep"].stop_timer();

      if (step % save_interval == 0 or step == Nt) {
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
