#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "program_params.hpp"
#include "particle_params.hpp"
#include "vec3.hpp"

#include <array>
#include <bitset>
#include <print>

namespace tf::particles {
  struct Particle {
    vec3<compute_t> location; // todo: use normalized locations here
    vec3<compute_t> old_location;
    vec3<compute_t> velocity;
    compute_t weight;
    double gamma;
  }; // end struct Particle

  struct alignas(alignof(Particle)) ParticleChunk {
    static constexpr std::size_t n_particles = 32;
    // (42 * 48 bytes + 3 * 4 bytes + 8 bytes + 4 bytes) = 2048 exactly
    //  ^ particles     ^ indices     ^ bitset  ^ padding

    // explicit ParticleChunk(const std::uint32_t i_, const std::uint32_t j_, const std::uint32_t k_) : idxs{i_, j_, k_} {}

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

    [[nodiscard]] std::size_t num_particles() const { return active.count(); }

    std::array<Particle, n_particles> particles{};
    std::bitset<n_particles> active{};
  }; // end struct ParticleChunk


  // struct ParticleChunkSoA {
  //   static constexpr std::size_t n_particles = 32;
  //   // (42 * 48 bytes + 3 * 4 bytes + 8 bytes + 4 bytes) = 2048 exactly
  //   //  ^ particles     ^ indices     ^ bitset  ^ padding
  //
  //   template<typename T> using array_t = std::array<T, n_particles>;
  //
  //   bool add_particle(const vec3<compute_t>& new_loc, const vec3<compute_t>& old_loc, const vec3<compute_t>& new_vel, const compute_t new_weight, const double new_gamma) {
  //     // todo: would disabling this check effect performance? Just loop over and then default false.
  //     if (active.all()) { return false; } // chunk is full
  //
  //     for (std::size_t i = 0; i < n_particles; ++i) {
  //       if (!active.test(i)) {
  //         location[i]     = new_loc;
  //         old_location[i] = old_loc;
  //         velocity[i]     = new_vel;
  //         weight[i]       = new_weight;
  //         gamma[i]        = new_gamma;
  //         active.set(i, true);
  //         return true;
  //       }
  //     }
  //     return false;
  //   }
  //
  //   void remove_particle(const std::size_t pid) {
  //     active.set(pid, false);
  //   }
  //
  //   [[nodiscard]] std::size_t num_particles() const { return active.count(); }
  //
  //   array_t<vec3<compute_t>> location;
  //   array_t<vec3<compute_t>> old_location;
  //   array_t<vec3<compute_t>> velocity;
  //   array_t<compute_t> weight;
  //   array_t<double> gamma;
  //   std::bitset<n_particles> active{};
  // }; // end struct Particle

  struct ParticleGroup {
    struct CellData {
      std::vector<ParticleChunk> chunks;
      std::array<std::uint32_t, 3> idxs{};
    };

    ParticleGroup() = delete;

    ParticleGroup(std::string name_, const compute_t mass_, const compute_t charge_, const std::size_t z_)
    : name(std::move(name_)),
      num_particles(0zu),
      atomic_number(z_),
      mass(mass_),
      charge(charge_),
      qdt_over_2m(static_cast<compute_t>(0.5 * static_cast<double>(charge) * static_cast<double>(dt) / static_cast<double>(mass))),
      cells{Ncx * Ncy * Ncz}
    {
      for (std::uint32_t i = 0; i < Ncx; i++) {
        for (std::uint32_t j = 0; j < Ncy; j++) {
          for (std::uint32_t k = 0; k < Ncz; k++) {
            const auto index = get_cid(i, j, k);
            cells[index] = {{}, {i, j, k}};
          } // end for(k)
        } // end for(j)
      } // end for(i)
    }

    void add_particle(const Particle& p, const std::size_t i, const std::size_t j, const std::size_t k) {
      const auto cid = get_cid(i, j, k);
      num_particles++;
      for (auto& chunk : cells[cid].chunks) {
        if (chunk.add_particle(p)) {
          return;
        } // particle added to chunk
      }
      // no available spots in current chunks
      cells[cid].chunks.emplace_back(); // add new chunk
      cells[cid].chunks.back().add_particle(p); // add particle to new chunk
    }

    void remove_particle(ParticleChunk& chunk, const std::size_t pid) {
      chunk.remove_particle(pid);
      num_particles--;
    }

    std::string name;
    std::size_t num_particles;
    std::size_t atomic_number;
    compute_t mass;
    compute_t charge;
    compute_t qdt_over_2m;

    std::vector<CellData> cells;
  }; // end struct ParticleGroup

  struct ParticleInitializer {
    static ParticleGroup initializeFromFile(const std::string&, compute_t, compute_t, std::size_t , const std::string&);
  };
} // end namespace tf::particles

#endif //PARTICLE_HPP
