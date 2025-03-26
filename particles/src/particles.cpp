#include "particles.hpp"

#include "program_params.hpp"

#include <print>

namespace tf::particles {
  void ParticleInitializer::initializeFromFile(ParticleGroup& g, const std::string& filename) {
    std::fstream file(filename, std::ios::in);

    if (!file.is_open()) {
      throw std::runtime_error("Particle initialization from file failed: " + filename);
    }

    vec3<compute_t> location{};
    vec3<compute_t> velocity{};
    float weight = 0.0;

    std::println("Opened particle file: {}", filename);
    std::string line;
    while (getline(file, line)) {
      std::istringstream buffer(line);
      buffer >> location >> velocity >> weight;// >> uid;

      const auto i = static_cast<std::size_t>(std::abs(location[0] / dx));
      const auto j = static_cast<std::size_t>(std::abs(location[1] / dy));
      const auto k = static_cast<std::size_t>(std::abs(location[2] / dz));

      // compute Lorentz factor and relativistic momentum
      double gamma = 1.0 / std::sqrt(1.0 - velocity.length_squared() / constants::c_sqr);
      auto momentum = static_cast<compute_t>(gamma) * g.mass * velocity;

      // add particle to group
      g.add_particle(Particle{location, location, momentum, weight, gamma}, i, j, k);
    }
    file.close();
  } // end initializeFromFile

}