#include "program_params.hpp"
#include "em_solver.hpp"
#include "metrics.hpp"
#include "timers.hpp"
#include "particles/particles.hpp"
#include "particles/pusher.hpp"
#include "particles/current_deposition.hpp"
#include "particles/collisions.hpp"

#include "barkeep.h"

using namespace tf;
using namespace tf::electromagnetics;
using namespace tf::particles;
using namespace tf::collisions;
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

   std::vector<ParticleGroup> particle_groups;
   for (const auto& name: particle_data) {
      ParticleInitializer::initializeFromFile(std::string{sim_path} + std::string{name}, particle_groups);
   }

   std::vector<Collisions> collisions;
   for (const auto& col : collision_params) {
      const auto& [name1, name2, clog, mult, step, self_scatter] = col;
      assert(!name1.empty() and !name2.empty()); // in case array has default initialized spots
      const auto g1 = std::ranges::find(particle_groups, name1, &ParticleGroup::name);
      const auto g2 = std::ranges::find(particle_groups, name2, &ParticleGroup::name);
      collisions.emplace_back(*g1, *g2, clog, mult, step, self_scatter);
   }

   emsolver_t emsolver(Nx, Ny, Nz);
   if constexpr (laser_enabled) {
      add_gaussianbeam(emsolver);
   }

   if constexpr (rmf_enabled) {
      add_rmf_antennas(emsolver, rmf_params);
   }

   emsolver.particle_correction();

   Metrics metrics(std::string{sim_path} + "/data/" + std::string{sim_name});
   if constexpr (em_enabled) {
      metrics.add_em_metrics(emsolver);
   }

   for (auto& g : particle_groups) {
      BorisPush::backstep_velocity(g, emsolver.emdata);
      metrics.add_particle_metric(g);
   }

   auto time = 0.0;
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
   metrics.write(step, time);
   timers["IO"].stop_timer();

   progress_bar->show();
   for (step = 1; step <= Nt; step++, time += dt) {
      // Electromagnetics
      timers["EM"].start_timer();
      emsolver.advance(time);
      timers["EM"].stop_timer();

      // Particle Push
      timers["Push"].start_timer();
      for (auto& g : particle_groups) {
         g.reset_positions();
         BorisPush::advance(g, emsolver.emdata, step);
      }
      timers["Push"].stop_timer();

      // Particle Sort (for timing)
      if (step % ParticleGroup::SORT_INTERVAL == 0) {
         timers["Sort"].start_timer();
         for (auto& g : particle_groups) {
            g.sort_particles();
         }
         timers["Sort"].stop_timer();
      }

      // Current Deposition
      timers["Jdep"].start_timer();
      for (auto& g : particle_groups) {
         CurrentDeposition::advance(g, emsolver.emdata);
      }
      timers["Jdep"].stop_timer();

      // Collisions
      timers["Collisions"].start_timer();
      for (auto& c : collisions) {
         c.advance(step);
      }
      timers["Collisions"].stop_timer();

      // Metrics output
      if (step % save_interval == 0 or step == Nt) {
         timers["IO"].start_timer();
         metrics.write(step, time);
         timers["IO"].stop_timer();
      }
   }
   progress_bar->done();
   timers["Main"].stop_timer();

   print_final_timers(timers);
   return 0;
}
