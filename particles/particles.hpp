#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "program_params.hpp"
#include "vec3.hpp"
#include "morton.hpp"
#include "constants.hpp"

// #include <gfx/timsort.hpp>
#include "dbg.h"

#include <vector>
#include <algorithm>
#include <print>
#include <fstream>
#include <sstream>

namespace tf::particles {
  struct Particle {
    vec3<compute_t> location; // todo: use normalized locations here
    vec3<compute_t> old_location;
    vec3<compute_t> velocity;
    compute_t weight;
    double gamma;
    bool disabled{false};
  }; // end struct Particle

  template<typename T = std::size_t>
  constexpr vec3<T> getCIDs(const vec3<compute_t>& loc) {
    return {
      static_cast<T>(std::floor(loc[0])),
      static_cast<T>(std::floor(loc[1])),
      static_cast<T>(std::floor(loc[2]))
    };
  }

  struct ParticleGroup {
    static constexpr std::size_t SORT_INTERVAL = 50;
    // static constexpr std::size_t DISABLED = morton_encode(Nx + 1, Ny + 1, Nz + 1);

    ParticleGroup() = delete;

    ParticleGroup(std::string name_,
                  const compute_t mass_,
                  const compute_t charge_,
                  const std::size_t z_)
    : name(std::move(name_)),
      atomic_number(z_),
      mass(mass_),
      charge(charge_),
      qdt_over_2m(calculate_qdt_over_2m())
    {}

    [[nodiscard]] compute_t calculate_qdt_over_2m() const {
      return static_cast<compute_t>(0.5 * static_cast<double>(charge) * static_cast<double>(dt) / static_cast<double>(mass));
    }

    [[nodiscard]] std::size_t num_particles() const { return particles.size(); }

    void reset_y_positions() {
      #pragma omp parallel for simd num_threads(nThreads)
      for (std::size_t pid = 0; pid < particles.size(); pid++) {
        particles[pid].location[1] = initial_y_position;
      }
    }

    void sort_particles() {
      // gfx::timsort(particles, {}, &Particle::code);
      std::ranges::sort(particles,
        [](const Particle& a, const Particle& b) {
          return morton_encode(getCIDs<std::size_t>(a.location)) < morton_encode(getCIDs<std::size_t>(b.location));
        }
      );
      std::erase_if(particles, [](const Particle& p) { return p.disabled; });
    }

    std::string name;
    std::size_t atomic_number;
    compute_t mass;
    compute_t charge;
    compute_t qdt_over_2m;
    compute_t initial_y_position{};
    std::vector<Particle> particles{};
  }; // end struct ParticleGroup

  struct ParticleInitializer {
    static ParticleGroup initializeFromFile(const std::string& name, const compute_t mass, const compute_t charge, const std::size_t z, const std::string& filename) {
      std::ifstream file(filename);

      if (!file.is_open()) {
        throw std::runtime_error("Particle initialization from file failed: " + filename);
      }

      ParticleGroup g(name, mass, charge, z);

      constexpr vec3 deltas{dx, dy, dz};
      constexpr vec3 mins{x_range[0], y_range[0], z_range[0]};

      vec3<double> location{};
      vec3<double> velocity{};
      float weight = 0.0;

      std::println("Loading particle file: {}... ", filename);
      std::string line;
      while (getline(file, line)) {
        std::istringstream buffer(line);
        buffer >> location >> velocity >> weight;

        location = (location - mins) / deltas;
        // compute Lorentz factor and relativistic momentum
        const auto gamma = 1.0 / std::sqrt(1.0 - velocity.length_squared() * constants::over_c_sqr<double>);

        // add particle to group
        g.particles.emplace_back(
          location.as_type<compute_t>(),
          location.as_type<compute_t>(),
          velocity.as_type<compute_t>() * gamma,
          weight,
          gamma,
          false
        );
      }
      file.close();
      g.sort_particles();
      g.initial_y_position = g.particles[0].location[1];
      return g;
    } // end initializeFromFile
  };
} // end namespace tf::particles

#endif //PARTICLE_HPP
