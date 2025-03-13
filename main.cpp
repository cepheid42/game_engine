#include "program_params.hpp"
#include "constants.hpp"
#include "em_solver.hpp"
#include "metrics.hpp"
#include "sources.hpp"
#include "array.hpp"

#include <print>

int main() {
  using array_t = tf::Array3D<double>;
  using ricker_t = tf::electromagnetics::RickerSource;
  using temporal_vec = std::vector<std::unique_ptr<tf::electromagnetics::TemporalSource>>;

  constexpr auto nx2 = Nx / 2;

  constexpr auto freq = tf::constants::c / (20.0 * dx);

  auto make_ricker = [&freq](temporal_vec& srcs) {
    srcs.push_back(std::make_unique<ricker_t>(freq));
  };

  auto make_srcvec = [&]() -> temporal_vec {
    temporal_vec result{};
    make_ricker(result);
    return result;
  };

  tf::electromagnetics::EMSolver emsolver(Nx, Ny, Nz, cfl, dt);

  emsolver.emdata.srcs.emplace_back(
    &emsolver.emdata.Ez,
    tf::electromagnetics::SpatialSource{
      make_srcvec(),
      1.0,
      {nx2, nx2 + 1, nx2, nx2 + 1, nx2, nx2 + 1}
    }
  );

  tf::metrics::Metrics metrics("/home/cepheid/TriForce/game_engine/data");

  metrics.addMetric(
    std::make_unique<tf::metrics::EMFieldsMetric>(
      std::unordered_map<std::string, array_t*> {
        {"Ez", &emsolver.emdata.Ez},
        // {"Ey", std::make_shared<array_t>(emsolver.emdata.Ey)}
      },
      metrics.adios.DeclareIO("EMFields")
    )
  );

  double t = 0.0;
  std::size_t step = 0zu;
  while (t < total_time) {
    emsolver.advance(t);

    if (step % 4 == 0) {
      const auto percent = 100.0 * t / total_time;
      std::println("Step {:4} Time: {:10.4e} Complete: {:4.1f}", step, t, percent);
      metrics.write(step);
    }

    t += dt;
    step++;
  }

  return 0;
}
