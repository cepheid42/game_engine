#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "program_params.hpp"
#include "vec3.hpp"

#include <array>
#include <bitset>
// #include <print>

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
    std::bitset<n_particles> moves{};
  }; // end struct ParticleChunk

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
      qdt_over_2m(calculate_qdt_over_2m()),
      cells{Ncx * Ncy * Ncz}
    {
      for (std::uint32_t i = 0; i < Ncx; i++) {
        for (std::uint32_t j = 0; j < Ncy; j++) {
          for (std::uint32_t k = 0; k < Ncz; k++) {
            cells[get_cid(i, j, k)] = {{}, {i, j, k}};
          } // end for(k)
        } // end for(j)
      } // end for(i)
    }

    [[nodiscard]] compute_t calculate_qdt_over_2m() const {
      return static_cast<compute_t>(0.5 * static_cast<double>(charge) * static_cast<double>(dt) / static_cast<double>(mass));
    }

    // void add_particle(const Particle& p, const std::size_t cid)
    void add_particle(const Particle& p, const std::size_t cid) {
;      num_particles++;
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

    void reset_y_positions() {
#pragma omp parallel for num_threads(nThreads)
      for (auto& [chunks, idxs] : cells) {
        if (chunks.empty()) { continue; }
        for (auto& chunk : chunks) {
          for (std::size_t i = 0; i < ParticleChunk::n_particles; ++i) {
            chunk[i].location[1] = initial_y_position;
          }
        }
      }
    }

    std::string name;
    std::size_t num_particles;
    std::size_t atomic_number;
    compute_t mass;
    compute_t charge;
    compute_t qdt_over_2m;
    compute_t initial_y_position{};

    std::vector<CellData> cells;
  }; // end struct ParticleGroup

  struct ParticleInitializer {
    static ParticleGroup initializeFromFile(const std::string&, compute_t, compute_t, std::size_t , const std::string&);
  };
} // end namespace tf::particles

#endif //PARTICLE_HPP
