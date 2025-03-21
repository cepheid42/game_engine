#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "program_params.hpp"
#include "vec3.hpp"
#include "constants.hpp"

#include <array>
#include <bitset>
#include <type_traits>
#include <string>
#include <fstream>

namespace tf::particles {
  struct Particle {
    vec3<compute_t> location;
    vec3<compute_t> old_location;
    vec3<compute_t> momentum;
    float weight;
    double gamma;
  }; // end struct Particle

  struct alignas(alignof(Particle)) ParticleChunk {
    using vec_t = vec3<compute_t>;

    static constexpr std::size_t padding = 40;
    static constexpr std::size_t n_particles = 84;
    // (84 * 48 bytes + 16 bytes + 40 bytes) = 4096 exactly
    //  ^ particles      ^ bitset   ^ padding

    explicit ParticleChunk(const std::size_t cid_) : cid(cid_) {}

    bool add_particle(const Particle& p) {
      // todo: would disabling this check effect performance? Just loop over and then default false.
      if (active.all()) { return false; } // chunk is full

      for (std::size_t i = 0; i < n_particles; ++i) {
        if (!active.test(i)) {
          particles[i] = p;
          active.set(i);
          return true;
        }
      }

      return false;
    }

    void remove_particle(const std::size_t pid) {
      active.set(pid, false);
    }

    Particle& operator[](const size_t pid) { return particles[pid]; }
    const Particle& operator[](const size_t pid) const { return particles[pid]; }

    [[nodiscard]] std::size_t num_active() const { return active.count(); }

    std::array<Particle, n_particles> particles{};
    std::bitset<n_particles> active{};
    std::uint32_t cid{};
    std::aligned_storage_t<40, alignof(Particle)> pad;
  }; // end struct ParticleChunk

  struct ParticleGroup {
    ParticleGroup() = delete;

    ParticleGroup(std::string name_, const compute_t mass_, const compute_t charge_, const std::size_t z_)
    : name(std::move(name_)),
      atomic_number(z_),
      mass(mass_),
      charge(charge_),
      inv_mass(1.0f / mass),
      cm_ratio(charge / mass),
      inv_cm_ratio(1.0f / cm_ratio),
      qdt_over_2m(0.5f * charge * dt * constants::q_e / mass),
      update_interval()
    {}

    [[nodiscard]] bool update_this_step(const compute_t time) {
      return update_interval.update_this_step(time);
    }

    void add_particle(Particle&& p, std::size_t cid) {
      // todo: this would be where setting octree bits would happen
      for (auto& chunk : cells[cid]) {
        if (chunk.add_particle(p)) {
          return;
        }
      }
      cells[cid].emplace_back(cid);
      cells[cid].back().add_particle(p);
    }


    std::string name;
    std::size_t atomic_number;
    compute_t mass;
    compute_t charge;
    compute_t inv_mass;
    compute_t cm_ratio;
    compute_t inv_cm_ratio;
    compute_t qdt_over_2m;

    Interval update_interval{};
    std::vector<std::vector<ParticleChunk>> cells{};
  }; // end struct ParticleGroup

  struct ParticleInitializer {
    static void initializeFromFile(ParticleGroup& g, const std::string& filename) {
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

        //std::cout << uid << group.mass << group.charge << std::endl;

        // compute Lorentz factor and relativistic momentum
        double gamma = 1.0 / std::sqrt(1.0 - velocity.length_squared() / constants::c_sqr);
        auto momentum = gamma * g.mass * velocity;

        // add particle to group
        g.add_particle({location, location, momentum, gamma, weight});
      }
      file.close();
    } // end initializeFromFile
  }; // end struct Particle Initializer
} // end namespace tf::particles

#endif //PARTICLE_HPP
