#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "program_params.hpp"
#include "vec3.hpp"
#include "morton.hpp"
#include "constants.hpp"

#include <gfx/timsort.hpp>

#include <vector>
#include <algorithm>
#include <print>
#include <fstream>

namespace tf::particles {
  struct Particle {
    vec3<compute_t> location; // todo: use normalized locations here
    vec3<compute_t> old_location;
    vec3<compute_t> velocity;
    compute_t weight;
    double gamma;
    std::size_t code;
  }; // end struct Particle


  struct ParticleGroup {
    static constexpr std::size_t SORT_INTERVAL = 50;
    static constexpr std::size_t DISABLED = morton_encode(Nx + 1, Ny + 1, Nz + 1);

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
      // #pragma omp parallel for simd num_threads(nThreads)
      for (std::size_t pid = 0; pid < particles.size(); pid++) {
        particles[pid].location[1] = initial_y_position;
      }
    }

    void sort_particles() {
      // gfx::timsort(particles, {}, &Particle::code);
      std::ranges::sort(particles, {}, &Particle::code);
      // std::ranges::stable_sort(particles, {}, &Particle::code);
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

      vec3 deltas{dx, dy, dz};
      vec3<double> location{};
      vec3<double> velocity{};
      float weight = 0.0;
      double y_init = 0.0;

      std::println("Loading particle file: {}... ", filename);
      std::string line;
      while (getline(file, line)) {
        std::istringstream buffer(line);
        buffer >> location >> velocity >> weight;

        const auto ix = static_cast<std::size_t>(std::abs((location[0] - x_range[0]) / dx));
        const auto iy = static_cast<std::size_t>(std::abs((location[1] - y_range[0]) / dy));
        const auto iz = static_cast<std::size_t>(std::abs((location[2] - z_range[0]) / dz));

        std::println("{}, {}, {}", ix, iy, iz);

        // for (std::size_t i = 0; i < 3; ++i) {
        //   const auto sx = location[i] / deltas[i];
        //   location[i] = sx - std::floor(sx);
        // }

        y_init = location[1];
        // compute Lorentz factor and relativistic momentum
        const auto gamma = 1.0 / std::sqrt(1.0 - velocity.length_squared() * constants::over_c_sqr<double>);

        // add particle to group
        g.particles.emplace_back(
          location.as_type<compute_t>(),
          location.as_type<compute_t>(),
          velocity.as_type<compute_t>(),
          weight,
          gamma,
          morton_encode(ix, iy, iz)
        );
      }
      file.close();
      g.sort_particles();
      g.initial_y_position = static_cast<compute_t>(y_init);
      return g;
    } // end initializeFromFile
  };
} // end namespace tf::particles

#endif //PARTICLE_HPP
