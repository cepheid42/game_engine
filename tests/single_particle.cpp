#include "program_params.hpp"
#include "constants.hpp"
#include "em_solver.hpp"
#include "metrics.hpp"
#include "array.hpp"
#include "timers.hpp"
#include "particles.hpp"
#include "pusher.hpp"

#include "barkeep.h"

#include <chrono>
#include <unordered_map>

using namespace tf;
using namespace tf::electromagnetics;
using namespace tf::particles;
using namespace tf::metrics;
namespace bk = barkeep;

Metrics create_metrics(const std::string& dir, EMSolver& em, const ParticleGroup& g1) {
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
   return metrics;
}


int main() {
   auto timers = utilities::create_timers();

   timers["Main"].start_timer();
   constexpr auto m_e = constants::m_e<compute_t>;
   constexpr auto q_e = constants::q_e<compute_t>;
   ParticleGroup g1("electrons", m_e, -q_e, 0);
   g1.particles.emplace_back(
      vec3{1.5_fp, 0.5_fp, 5.5_fp}, // location
      vec3{1.5_fp, 0.5_fp, 5.5_fp}, // old location
      vec3{0.0_fp, 0.0_fp, 0.0_fp}, // velocity
      1.0_fp,                       // weight
      1.0_fp,                       // gamma
      false                         // disabled
   );

   g1.initial_y_position = 0.5_fp;

   EMSolver emsolver(Nx, Ny, Nz, cfl, dt);
   emsolver.emdata.Ex_app.fill(-1000.0_fp);

   constexpr BorisPush particle_push{};
   BorisPush::backstep_velocity(g1, emsolver.emdata);

   const auto metrics = create_metrics(
      "/home/cepheid/TriForce/game_engine/data/single_particle_test",
      emsolver,
      g1
   );

   std::size_t step = 0zu;
   const auto progress_bar =
      bk::ProgressBar(
         &step, {
            .total = Nt,
            .message = "Step",
            .speed = 0.,
            .speed_unit = "steps/s",
            .interval = 1.,
            // .no_tty = true,
            .show = false
         });

   timers["IO"].start_timer();
   metrics.write(step);
   timers["IO"].stop_timer();

   progress_bar->show();
   for (std::size_t n = 0; n < Nt; n++) {
      timers["EM"].start_timer();
      emsolver.advance(static_cast<compute_t>(n) * dt);
      timers["EM"].stop_timer();

      timers["Push"].start_timer();
      particle_push(g1, emsolver.emdata, step);
      timers["Push"].stop_timer();

      g1.reset_y_positions();

      step++;
      if (step % save_interval == 0) {
         timers["IO"].start_timer();
         metrics.write(step);
         timers["IO"].stop_timer();
      }
   }
   timers["Main"].stop_timer();
   progress_bar->done();

   print_final_timers(timers);
   return 0;
}