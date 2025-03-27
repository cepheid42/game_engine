#include "particles.hpp"

#include "program_params.hpp"

#include <print>

namespace tf::particles {
  void ParticleInitializer::initializeFromFile(ParticleGroup& g, const std::string& filename) {
    std::fstream file(filename, std::ios::in);

    if (!file.is_open()) {
      throw std::runtime_error("Particle initialization from file failed: " + filename);
    }

    vec3<double> deltas{dx, dy, dz};
    vec3<double> location{};
    vec3<double> velocity{};
    float weight = 0.0;

    std::println("Opened particle file: {}", filename);
    std::string line;
    while (getline(file, line)) {
      std::istringstream buffer(line);
      buffer >> location >> velocity >> weight;// >> uid;

      std::array<std::size_t, 3> index{
        static_cast<std::size_t>(std::abs(location[0] / dx)),
        static_cast<std::size_t>(std::abs(location[1] / dy)),
        static_cast<std::size_t>(std::abs(location[2] / dz))
      };

      for (std::size_t i = 0; i < 3; ++i) {
        const auto sx = location[i] / deltas[i];
        location[i] = sx - std::floor(sx);
      }

      // compute Lorentz factor and relativistic momentum
      double gamma = 1.0 / std::sqrt(1.0 - velocity.length_squared() / constants::c_sqr);

      // add particle to group
      g.add_particle(
        Particle{
          location.to_float(),
          location.to_float(),
          velocity.to_float(),
          weight,
          gamma,
          true
        },
        index);
    }
    file.close();
  } // end initializeFromFile

}