// #include "program_params.hpp"
// #include "constants.hpp"
// #include "em_solver.hpp"
// #include "metrics.hpp"
// #include "sources.hpp"
// #include "array.hpp"
// #include "timers.hpp"
#include "particle.hpp"
// #include "pusher.hpp"
//

#include "octree.hpp"

// #include <print>
// #include <bitset>
// #include <iostream>
// #include <chrono>

int main() {
  tf::particles::ParticleGroup p("test", 0.0, 0.0, 1);

  p.cells[250].push_back(tf::particles::Particle{{}, {}, {}, 250.0, 0.0});
  p.cells[251].push_back(tf::particles::Particle{{}, {}, {}, 251.0, 0.0});
  p.cells[252].push_back(tf::particles::Particle{{}, {}, {}, 252.0, 0.0});
  p.cells[253].push_back(tf::particles::Particle{{}, {}, {}, 253.0, 0.0});
  p.cells[254].push_back(tf::particles::Particle{{}, {}, {}, 254.0, 0.0});

  p.update_tree();

  tf::particles::visit_octree(p.tree);

  // for (const auto& c : p.cells) {
  //   for (const auto& ps : c) {
  //     std::println("{}", ps.weight);
  //   }
  // }

  // for (std::size_t i = 0; i < 8; i++) {
  //   for (std::size_t j = 0; j < 8; j++) {
  //     for (std::size_t k = 0; k < 8; k++) {
  //       auto idx = tf::particles::morton_encode(i, j, k);
  //
  //     }
  //   }
  // }

  // std::println("Particle size: {} bytes", sizeof(tf::particles::Particle));
  // std::println("bitset<85> size: {} bytes", sizeof(std::bitset<85>));
  // std::println("Chunk size: {} bytes", sizeof(tf::particles::ParticleChunk), sizeof(tf::particles::Particle), sizeof(std::bitset<64>));
  //

  // using array_t = tf::Array3D<compute_t>;
  // using ricker_t = tf::electromagnetics::RickerSource;
  // using continuous_t = tf::electromagnetics::ContinuousSource;
  // using temporal_vec = std::vector<std::unique_ptr<tf::electromagnetics::TemporalSource>>;
  //
  // constexpr auto nx2 = Nx / 2;
  // constexpr auto freq = tf::constants::c / (20.0f * dx);
  //
  // // auto make_ricker = [&freq](temporal_vec& srcs) {
  // //   srcs.push_back(std::make_unique<ricker_t>(freq));
  // // };
  //
  // auto make_continuous = [&freq](temporal_vec& srcs) {
  //   srcs.push_back(std::make_unique<continuous_t>(2.0 * tf::constants::pi * freq, 0.0f, 0.0f, 1.0e30f, dx));
  // };
  //
  // auto make_srcvec = [&]() -> temporal_vec {
  //   temporal_vec result{};
  //   // make_ricker(result);
  //   make_continuous(result);
  //   return result;
  // };
  //
  // tf::electromagnetics::EMSolver emsolver(Nx, Ny, Nz, cfl, dt);
  //
  // emsolver.emdata.srcs.emplace_back(
  //   &emsolver.emdata.Ez,
  //   tf::electromagnetics::SpatialSource{
  //     make_srcvec(),
  //     100.0f,
  //     {nx2, nx2 + 1, nx2, nx2 + 1, nx2, nx2 + 1}
  //   }
  // );
  //
  // tf::metrics::Metrics metrics("/home/cepheid/TriForce/game_engine/data/single_continuous_test");
  //
  // metrics.addMetric(
  //   std::make_unique<tf::metrics::EMFieldsMetric>(
  //     std::unordered_map<std::string, array_t*> {
  //       {"Ex", &emsolver.emdata.Ex},
  //       {"Ey", &emsolver.emdata.Ey},
  //       {"Ez", &emsolver.emdata.Ez},
  //       {"Hx", &emsolver.emdata.Hx},
  //       {"Hy", &emsolver.emdata.Hy},
  //       {"Hz", &emsolver.emdata.Hz},
  //     },
  //     metrics.adios.DeclareIO("EMFields")
  //   )
  // );
  //
  // Timer timer{};
  //
  // // std::println("dt = {}", dt);
  //
  // compute_t t = 0.0f;
  // std::size_t step = 0zu;
  //
  // timer.start_timer();
  // while (t <= total_time) {
  //   emsolver.advance(t);
  //
  //   if (step % save_interval == 0) {
  //     const auto percent = 100.0f * t / total_time;
  //     std::println("Step {:4} Time: {:7.1e} Complete: {:4.1f}%", step, t, percent);
  //     metrics.write(step);
  //   }
  //
  //   t += dt;
  //   step++;
  // }
  // const auto percent = 100.0f * t / total_time;
  // std::println("Step {:4} Time: {:7.1e} Complete: {:4.1f}%", step, t, percent);
  // metrics.write(step);
  //
  // timer.stop_timer();
  //
  // const auto time = std::chrono::hh_mm_ss(timer.elapsed);
  // std::println("Total Time: {}", time);
  return 0;
}
