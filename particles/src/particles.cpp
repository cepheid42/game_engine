#include "particles.hpp"

#include "program_params.hpp"

#include <print>
#include <fstream>

namespace tf::particles {
  ParticleGroup ParticleInitializer::initializeFromFile(const std::string& name, const compute_t mass, const compute_t charge, const std::size_t z, const std::string& filename) {
    std::ifstream file(filename);

    if (!file.is_open()) {
      throw std::runtime_error("Particle initialization from file failed: " + filename);
    }

    ParticleGroup g(name, mass, charge, z);

    vec3<double> deltas{dx, dy, dz};
    vec3<double> location{};
    vec3<double> velocity{};
    float weight = 0.0;

    std::println("Opened particle file: {}", filename);
    std::string line;
    while (getline(file, line)) {
      std::istringstream buffer(line);
      buffer >> location >> velocity >> weight;// >> uid;

      // std::println("{}, {}, {}", location[0], location[1], location[2]);

      std::array<std::size_t, 3> index{
        static_cast<std::size_t>(std::abs((location[0] - x_range[0]) / dx)),
        static_cast<std::size_t>(std::abs((location[1] - y_range[0]) / dy)),
        static_cast<std::size_t>(std::abs((location[2] - z_range[0]) / dz))
      };

      // std::println("{}, {}, {} = {}", index[0], index[1], index[2], morton_encode(index[0], index[1], index[2]));

      for (std::size_t i = 0; i < 3; ++i) {
        const auto sx = location[i] / deltas[i];
        location[i] = sx - std::floor(sx);
      }

      // compute Lorentz factor and relativistic momentum
      const auto gamma = 1.0 / std::sqrt(1.0 - velocity.length_squared() * constants::over_c_sqr);

      // add particle to group
      g.add_particle(
        Particle{
          location.as_type<compute_t>(),
          location.as_type<compute_t>(),
          (velocity * gamma).as_type<compute_t>(),
          weight,
          gamma,
          true
        },
        index);
    }
    file.close();

    g.update_tree();
    return g;
  } // end initializeFromFile

}