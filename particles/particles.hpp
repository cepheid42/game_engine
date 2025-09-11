#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "program_params.hpp"
#include "vec3.hpp"
#include "morton.hpp"
#include "constants.hpp"

// #include <gfx/timsort.hpp>

#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>

namespace tf::particles {
struct Particle {
   vec3<double> location;
   vec3<double> old_location;
   vec3<double> velocity;
   double weight;
   double gamma;
   bool disabled{false};
}; // end struct Particle

template <typename T = std::size_t>
constexpr vec3<T> getCellIndices(const vec3<double>& loc) {
   return {
      static_cast<T>(std::floor(loc[0])),
      static_cast<T>(std::floor(loc[1])),
      static_cast<T>(std::floor(loc[2]))
   };
}

constexpr std::size_t getCellIndex(const vec3<double>& loc) {
   const auto x = static_cast<std::size_t>(std::floor(loc[0]));
   const auto y = static_cast<std::size_t>(std::floor(loc[1]));
   const auto z = static_cast<std::size_t>(std::floor(loc[2]));
   return z + (Ncz * y) + (Ncy * Ncz * x);
}

constexpr auto calculateGamma(const auto& v) {
   // Calculates gamma using regular velocity
   return 1.0 / std::sqrt(1.0 - v.length_squared() * constants::over_c_sqr<double>);
}

constexpr auto calculateGammaV(const auto& v) {
   // Calculates gamma using gamma*v (e.g. relativistic momentum but with mass terms canceled)
   return std::sqrt(1.0 + v.length_squared() * constants::over_c_sqr<double>);
}

struct ParticleGroup {
   static constexpr std::size_t SORT_INTERVAL = 50;
   std::string name;
   double mass;
   double charge;
   double qdt_over_2m;
   // double initial_y_position{};
   std::vector<Particle> particles{};

   ParticleGroup() = delete;

   ParticleGroup(std::string name_, const double mass_, const double charge_)
   : name(std::move(name_)),
     mass(mass_),
     charge(charge_),
     qdt_over_2m(0.5 * charge * dt / mass)
   {}

   [[nodiscard]] std::size_t num_particles() const { return particles.size(); }

   // void reset_y_positions() {
   //    #pragma omp parallel for simd num_threads(nThreads)
   //    for (std::size_t pid = 0; pid < particles.size(); pid++) {
   //       particles[pid].location[1] = initial_y_position;
   //    }
   // }

   void sort_particles() {
      std::erase_if(particles, [](const Particle& p) { return p.disabled; });
      // gfx::timsort(
      //    particles,
      //    [](const Particle& a, const Particle& b) {
      //       return morton_encode(getCellIndices<std::size_t>(a.location)) < morton_encode(getCellIndices<std::size_t>(b.location));
      //    }
      // );
      std::ranges::sort(
         particles,
         [](const Particle& a, const Particle& b) {
            return morton_encode(getCellIndices<std::size_t>(a.location)) < morton_encode(getCellIndices<std::size_t>(b.location));
         }
      );
   }
}; // end struct ParticleGroup

struct ParticleInitializer {
   static ParticleGroup initializeFromFile(const std::string& filename) {
      std::ifstream file(filename);

      if (!file.is_open()) {
         throw std::runtime_error("Particle initialization from file failed: " + filename);
      }

      constexpr vec3 deltas{dx, dy, dz};
      constexpr vec3 mins{x_range[0], y_range[0], z_range[0]};

      std::string comment{};
      std::string name{};
      double mass{};
      double charge{};
      vec3<double> location{};
      vec3<double> velocity{};
      double weight = 0.0;

      std::println("Loading particle file: {}... ", filename);

      std::string line;
      std::getline(file, line);
      std::istringstream buff(line);

      buff >> comment >> name >> mass >> charge;

      ParticleGroup g(name, mass, charge);
      while (std::getline(file, line)) {
         std::istringstream buffer(line);
         buffer >> location >> velocity >> weight;

         location = (location - mins) / deltas;

         // compute Lorentz factor and relativistic momentum
         const auto gamma = calculateGamma(velocity);

         // add particle to group
         g.particles.emplace_back(
            location,
            location,
            velocity,
            weight,
            gamma,
            false
         );
      }
      file.close();
      if (g.particles.empty()) {
         throw std::runtime_error("Particle initialization failed: Particles vector is empty.");
      }
      g.sort_particles();
      // g.initial_y_position = g.particles.empty() ? 0.0 : g.particles[0].location[1];
      return g;
   } // end initializeFromFile
}; // end struct ParticleInitializer
} // end namespace tf::particles

#endif //PARTICLE_HPP
