#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "program_params.hpp"
#include "particle_params.hpp"
#include "vec3.hpp"
#include "constants.hpp"
#include "octree.hpp"

#include <fstream>
#include <print>

// template<size_t N>
// struct std::formatter<std::bitset<N>> : std::formatter<std::string> {
//   auto format(std::bitset<N> b, format_context& ctx) const {
//     return formatter<std::string>::format(b.to_string(), N);
//   }
// };

namespace tf::particles {
  struct Particle {
    vec3<compute_t> location; // todo: use normalized locations here
    vec3<compute_t> old_location;
    vec3<compute_t> velocity;
    compute_t weight;
    // double gamma;
    bool active;
  }; // end struct Particle

  struct ParticleCell {
    std::vector<Particle> particles;
    std::size_t cid;
  };

  struct ParticleGroup {
    using p_vec = std::vector<Particle>;

    ParticleGroup() = delete;

    ParticleGroup(std::string name_, const compute_t mass_, const compute_t charge_, const std::size_t z_)
    : name(std::move(name_)),
      num_particles(0zu),
      atomic_number(z_),
      mass(mass_),
      charge(charge_),
      // inv_mass(1.0f / mass),
      // cm_ratio(charge / mass),
      // inv_cm_ratio(1.0f / cm_ratio),
      qdt_over_2m(static_cast<compute_t>(0.5 * static_cast<double>(charge) * static_cast<double>(dt) / static_cast<double>(mass))),
      cells{(Nx - 1) * (Ny - 1) * (Nz - 1)},
      tree(create_particle_octree(cells))
    {
      for (std::size_t i = 0; i < Nx - 1; i++) {
        for (std::size_t j = 0; j < Ny - 1; j++) {
          for (std::size_t k = 0; k < Nz - 1; k++) {
            const auto code = morton_encode(i, j, k);
            cells[code] = {{}, code};
          }
        }
      }

    }

    void add_particle(Particle&& p, const std::size_t cid) {
      cells[cid].particles.push_back(p);
      num_particles++;
    }

    void add_particle(Particle&& p, const std::array<std::size_t, 3>& idx) {
      add_particle(std::forward<Particle>(p), morton_encode(idx[0], idx[1], idx[2]));
    }

    static bool update_tree_nodes(auto& node) {
      for (std::size_t i = 0; i < 8; i++) {
        if (node.is_leaf) {
          // check for particles in each cell and set the appropriate active bit
          node.active.set(i, !node.cells[i]->particles.empty());
        } else {
          // recurse
          const auto has_particles = update_tree_nodes(node.children[i]);
          node.active.set(i, has_particles);
        }
      }

      return node.active.any();
    }

    void update_tree() {
      update_tree_nodes(tree);
    }

    std::string name;
    std::size_t num_particles;
    std::size_t atomic_number;
    compute_t mass;
    compute_t charge;
    // compute_t inv_mass;
    // compute_t cm_ratio;
    // compute_t inv_cm_ratio;
    compute_t qdt_over_2m;

    std::vector<ParticleCell> cells;
    Octree<ParticleCell> tree;
  }; // end struct ParticleGroup

  struct ParticleInitializer {
    static ParticleGroup initializeFromFile(const std::string&, compute_t, compute_t, std::size_t , const std::string&);
  };
} // end namespace tf::particles

#endif //PARTICLE_HPP
