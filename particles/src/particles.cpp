#include "particles.hpp"

#include "program_params.hpp"
#include "constants.hpp"

#include <print>
#include <fstream>

namespace tf::particles {
  ParticleGroup ParticleInitializer::initializeFromFile(const std::string& name, const compute_t mass, const compute_t charge, const std::size_t z, const std::string& filename) {
    std::ifstream file(filename);

    if (!file.is_open()) {
      throw std::runtime_error("Particle initialization from file failed: " + filename);
    }

    ParticleGroup g(name, mass, charge, z);

    vec3 deltas{dx, dy, dz};
    vec3<double> location{};
    vec3<double> velocity{};
    float weight = 0.0;
    double y_init = 0.0;

    std::println("Opened particle file: {}", filename);
    std::string line;
    while (getline(file, line)) {
      std::istringstream buffer(line);
      buffer >> location >> velocity >> weight;

      const auto ix = static_cast<std::size_t>(std::abs((location[0] - x_range[0]) / dx));
      const auto iy = static_cast<std::size_t>(std::abs((location[1] - y_range[0]) / dy));
      const auto iz = static_cast<std::size_t>(std::abs((location[2] - z_range[0]) / dz));

      for (std::size_t i = 0; i < 3; ++i) {
        const auto sx = location[i] / deltas[i];
        location[i] = sx - std::floor(sx);
      }

      y_init = location[1];
      // compute Lorentz factor and relativistic momentum
      const auto gamma = 1.0 / std::sqrt(1.0 - velocity.length_squared() * constants::over_c_sqr<double>);

      const auto p = Particle{
        location.as_type<compute_t>(),
        location.as_type<compute_t>(),
        velocity.as_type<compute_t>(),
        weight,
        gamma
      };
      // add particle to group
      g.add_particle(p, get_cid(ix, iy, iz));
      // g.add_particle(p, get_cid(ix, iy, iz));
    }
    file.close();
    g.initial_y_position = static_cast<compute_t>(y_init);
    return g;
  } // end initializeFromFile
} // end namespace tf::particles