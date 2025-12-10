//
// #include "mdspan.hpp"
//
// #include <vector>
// #include <print>
// #include <cassert>
//
//
// // Full 3d span over data
// using mdspan_t = std::mdspan<
//    std::size_t,
//    std::dextents<std::size_t, 3>
// >;
//
// // 3d span
// using subspan3d_t = std::mdspan<
//    std::size_t,
//    std::dextents<std::size_t, 3>,
//    std::layout_stride
// >;
//
// // 2D span
// using subspan2d_t = std::mdspan<
//    std::size_t,
//    std::dextents<std::size_t, 2>,
//    std::layout_stride
// >;
//
// // 1D span
// using subspan_t = std::mdspan<
//    std::size_t,
//    std::dextents<std::size_t, 1>,
//    std::layout_stride
// >;
//
// void fill(auto& md) {
//    for (std::size_t i = 0; i < md.extent(0); i++) {
//       for (std::size_t j = 0; j < md.extent(1); j++) {
//          for (std::size_t k = 0; k < md.extent(2); k++) {
//             md[i, j, k] = k  + (md.extent(2) * j) + (md.extent(1) * md.extent(2) * i);
//          }
//       }
//    }
// }
//
// void difference(auto& src, auto& dest) {
//    for (std::size_t i = 0; i < dest.extent(0); i++) {
//       for (std::size_t j = 0; j < dest.extent(1); j++) {
//          for (std::size_t k = 0; k < dest.extent(2); k++) {
//             dest[i, j, k] = src[i, j, k + 1] - src[i, j, k];
//          }
//       }
//    }
// }
//
// void print(auto& md) {
//    for (std::size_t i = 0; i < md.extent(0); i++) {
//       for (std::size_t j = 0; j < md.extent(1); j++) {
//          for (std::size_t k = 0; k < md.extent(2); k++) {
//             std::print("{:>3}, ", md[i, j, k]);
//          }
//          std::println();
//       }
//       std::println();
//    }
//    std::println();
// }
// constexpr auto nx = 10zu;
// constexpr auto ny = 10zu;
// constexpr auto nz = 10zu;
// constexpr auto size = nx * ny * nz;
// int main() {
//    // Data
//    std::vector<std::size_t> test(size);
//    std::vector<std::size_t> result_data(size);
//    // Using mdspan make indexing subspans easier
//    mdspan_t og{test.data(), std::extents{nx, ny, nz}};
//    // mdspan_t result{result_data.data(), std::extents{nx, ny, nz}};
//
//    // fill(og);
//    // print(og);
//
//    subspan3d_t slice{&og[0, 1, 1], {std::extents{nx, ny - 2, nz - 2}, std::array{ny * nz, nz, 1zu}}};
//
//    fill(slice);
//    print(og);
//
//    return 0;
// }


#include "program_params.hpp"
#include "em_solver.hpp"
#include "em_data.hpp"
#include "metrics.hpp"
#include "timers.hpp"
// #include "particles/particles.hpp"
// #include "particles/pusher.hpp"
// #include "particles/current_deposition.hpp"
// #include "particles/collisions.hpp"

#include "barkeep.h"

#include <print>

using namespace tf;
using namespace tf::electromagnetics;
// using namespace tf::particles;
// using namespace tf::collisions;
// using namespace tf::metrics;

namespace bk = barkeep;


int main() {
   auto timers = utilities::create_timers();
   timers["Main"].start_timer();

   // std::vector<ParticleGroup> particle_groups;
   // for (const auto& name: particle_data) {
   //    ParticleInitializer::initializeFromFile(std::string{sim_path} + std::string{name}, particle_groups);
   // }

   // std::vector<Collisions> collisions;
   // for (const auto& col : collision_params) {
   //    const auto& [name1, name2, clog, mult, step, self_scatter] = col;
   //    assert(!name1.empty() and !name2.empty()); // in case array has default initialized spots
   //    const auto g1 = std::ranges::find(particle_groups, name1, &ParticleGroup::name);
   //    const auto g2 = std::ranges::find(particle_groups, name2, &ParticleGroup::name);
   //    collisions.emplace_back(*g1, *g2, clog, mult, step, self_scatter);
   // }

   const auto data_path = std::string{sim_path} + "/data/" + std::string{sim_name};
   auto emdata = make_emdata();
   // add_gaussianbeam(emdata);

   // EMSolver::particle_correction(emdata);
   // for (auto& g : particle_groups) {
   //    BorisPush::backstep_velocity(g, emdata);
   // }

   timers["IO"].start_timer();
   size_t padding = 10;
   std::string count_padded = std::to_string(0);
   count_padded.insert(count_padded.begin(), padding - count_padded.length(), '0');
   count_padded += ".bp"; // default extension is .bp

   metrics::EMFieldsMetric::write(emdata, data_path, count_padded);
   // metrics::ParticleDumpMetric::write(particle_groups[0], data_path, count_padded);
   // metrics::ParticleDumpMetric::write(particle_groups[1], data_path, count_padded);
   timers["IO"].stop_timer();


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


   progress_bar->show();
   for (step = 1; step <= Nt; step++, time += dt) {
      // std::println("Step {}", step);

      // Electromagnetics
      timers["EM"].start_timer();
      EMSolver::advance(emdata, step);
      timers["EM"].stop_timer();

      // // Particle Push
      // timers["Push"].start_timer();
      // for (auto& g : particle_groups) {
      //    if (step % 50 == 0) { g.sort_particles(); }
      //    g.reset_positions();
      //    BorisPush::advance(g, emdata);
      // }
      // timers["Push"].stop_timer();

      // // Current Deposition
      // timers["Jdep"].start_timer();
      // for (auto& g : particle_groups) {
      //    CurrentDeposition::advance(g, emsolver.emdata);
      // }
      // timers["Jdep"].stop_timer();
      //
      // // Collisions
      // timers["Collisions"].start_timer();
      // for (auto& c : collisions) {
      //    c.advance(step);
      // }
      // timers["Collisions"].stop_timer();

      timers["IO"].start_timer();
      if (step % em_save_interval == 0) {
         padding = 10;
         count_padded = std::to_string(step);
         count_padded.insert(count_padded.begin(), padding - count_padded.length(), '0');
         count_padded += ".bp"; // default extension is .bp

         metrics::EMFieldsMetric::write(emdata, data_path, count_padded);
         // metrics::ParticleDumpMetric::write(particle_groups[0], data_path, count_padded);
         // metrics::ParticleDumpMetric::write(particle_groups[1], data_path, count_padded);
      }
      timers["IO"].stop_timer();
   }
   progress_bar->done();
   timers["Main"].stop_timer();

   print_final_timers(timers);
   return 0;
}
