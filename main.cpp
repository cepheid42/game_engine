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
#include "hillsvortex.hpp"

#include "barkeep.h"

// #include <print>

using namespace tf;
using namespace tf::electromagnetics;
using namespace tf::particles;
using namespace tf::metrics;

namespace bk = barkeep;

Metrics create_metrics(const std::string& dir, EMSolver& em, const ParticleGroup& g1){//, const ParticleGroup& g2) {
   Metrics metrics(dir);

   metrics.addMetric(
      std::make_unique<EMFieldsMetric>(
         std::unordered_map<std::string, Array3D<compute_t>*>{
            {"Ex", &em.emdata.Ex_total},
            {"Ey", &em.emdata.Ey_total},
            {"Ez", &em.emdata.Ez_total},
            {"Hx", &em.emdata.Bx_total},
            {"Hy", &em.emdata.By_total},
            {"Hz", &em.emdata.Bz_total},
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

   // metrics.addMetric(
   //    std::make_unique<ParticleDumpMetric>(
   //       &g2,
   //       metrics.adios.DeclareIO(g2.name + "_dump")
   //    )
   // );
   //
   // metrics.addMetric(
   //    std::make_unique<ParticleMetric>(
   //       &g2,
   //       metrics.adios.DeclareIO(g2.name + "_metrics"),
   //       Ncx, Ncy, Ncz
   //    )
   // );

   return metrics;
}


int main() {
   auto timers = utilities::create_timers();
   timers["Main"].start_timer();
   constexpr auto q_e = constants::q_e<compute_t>;
   constexpr auto m_e = constants::m_e<compute_t>;
   // constexpr auto m_p = constants::m_p<compute_t>;
   // constexpr auto electron_file = "/home/cepheid/TriForce/game_engine/data/electron_slab.dat";
   // constexpr auto      ion_file = "/home/cepheid/TriForce/game_engine/data/ion_slab.dat";
   //
   // auto g1 = ParticleInitializer::initializeFromFile("electrons", m_e, -q_e, 0, electron_file);
   // auto g2 = ParticleInitializer::initializeFromFile(     "ions", m_p, +q_e, 1,      ion_file);

   ParticleGroup g1("electrons", m_e, -q_e, 0);

   constexpr vec3 loc{30.5, 0.5, 50.5};
   constexpr vec3 vel{0.0, 0.0, 1.87e7};

   g1.particles.emplace_back(
      loc,
      loc,
      vel,
      1.0,
      calculateGamma(vel),
      false
   );

   EMSolver emsolver(Nx, Ny, Nz);

   constexpr auto B0 = 100000.0;
   constexpr auto rs = 0.05;
   constexpr auto kk = 3.0;
   HillsFieldGenerator::fillEMData(emsolver.emdata, B0, rs, kk);

   for (std::size_t i = 0; i < emsolver.emdata.By_app.nx(); i++) {
      for (std::size_t j = 0; j < emsolver.emdata.By_app.ny(); j++) {
         for (std::size_t k = 0; k < emsolver.emdata.By_app.nz(); k++) {
            emsolver.emdata.By_app(i, j, k) = static_cast<compute_t>(i) * 0.12 * dx / (x_range[1] - x_range[0]);
         }
      }
   }

   emsolver.particle_correction();
   // add_gaussianbeam(emsolver);

   BorisPush::backstep_velocity(g1, emsolver.emdata);
   // BorisPush::backstep_velocity(g2, emsolver.emdata);

   const auto metrics = create_metrics(
      "/home/cepheid/TriForce/game_engine/data/lsi_test",
      emsolver,
      g1//, g2
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
          // .no_tty = true,
          .show = false});

   timers["IO"].start_timer();
   metrics.write(step);
   timers["IO"].stop_timer();

   progress_bar->show();
   for (step = 1; step <= Nt; step++) {
      // if (step == 10) exit(0);

      timers["EM"].start_timer();
      emsolver.advance(static_cast<compute_t>(step) * dt);
      timers["EM"].stop_timer();

      timers["Push"].start_timer();
      g1.reset_y_positions();
      // g2.reset_y_positions();
      BorisPush::advance(g1, emsolver.emdata, step);
      // BorisPush::advance(g2, emsolver.emdata, step);
      timers["Push"].stop_timer();

      // timers["Jdep"].start_timer();
      // CurrentDeposition::advance(g1, emsolver.emdata);
      // CurrentDeposition::advance(g2, emsolver.emdata);
      // timers["Jdep"].stop_timer();

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
