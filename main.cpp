#include "program_params.hpp"
// #include "em_solver.hpp"
#include "metrics.hpp"
#include "timers.hpp"
#include "particles/particles.hpp"
// #include "particles/pusher.hpp"
// #include "particles/current_deposition.hpp"
#include "particles/collisions.hpp"

#include "barkeep.h"

#include <ranges>
#include <unordered_map>

using namespace tf;
// using namespace tf::electromagnetics;
using namespace tf::particles;
using namespace tf::collisions;
using namespace tf::metrics;

namespace bk = barkeep;

int main() {
   std::println("Beginning simulation: \"{}\"...", sim_name);

   auto timers = utilities::create_timers();
   timers["Main"].start_timer();

   std::unordered_map<std::string, ParticleGroup> particle_groups{};
   for (const auto& species: particle_spec) {
      particle_groups.insert({std::string{species.name}, ParticleGroup(species)});
   }

   std::vector<Collisions> collisions;
   for (const auto& col : collision_spec) {
      collisions.emplace_back(col, particle_groups);
   }

   // emsolver_t emsolver(Nx, Ny, Nz);
   // add_gaussianbeam(emsolver);


   Metrics metrics(std::string{sim_path} + "/data/" + std::string{sim_name});
   // if constexpr (em_enabled) {
   //    metrics.add_em_metrics(emsolver);
   // }

   if constexpr (push_enabled or coll_enabled) {
      // emsolver.particle_correction();
      for (auto& g : particle_groups | std::views::values) {
         // BorisPush::backstep_velocity(g, emsolver.emdata);
         metrics.add_particle_metric(g);
      }
   }

   auto time = 0.0;
   auto step = 0zu;
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
      // std::println("--------------- Step {} ---------------", step);

      // // Electromagnetics
      // timers["EM"].start_timer();
      // emsolver.advance(time);
      // timers["EM"].stop_timer();
      //
      // // Particle Push
      // timers["Push"].start_timer();
      // for (auto& g : particle_groups | std::views::values) {
      //    g.reset_positions();
      //    BorisPush::advance(g, emsolver.emdata, step);
      // }
      // timers["Push"].stop_timer();
      //
      // // Current Deposition
      // timers["Jdep"].start_timer();
      // for (auto& g : particle_groups) {
      //    CurrentDeposition::advance(g, emsolver.emdata);
      // }
      // timers["Jdep"].stop_timer();

      // Collisions
      timers["Collisions"].start_timer();
      for (auto& c : collisions) {
         c.advance(step);
      }
      timers["Collisions"].stop_timer();

      // Metrics output
      timers["IO"].start_timer();
      metrics.write(step, time);
      timers["IO"].stop_timer();
   }
   progress_bar->done();
   timers["Main"].stop_timer();

   print_final_timers(timers);
   return 0;
}
