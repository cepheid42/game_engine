#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "program_params.hpp"
#include "particle_params.hpp"
#include "vec3.hpp"
#include "constants.hpp"
#include "octree.hpp"

// #include <array>
// #include <bitset>
// #include <type_traits>
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
      // cells{Nx * Ny * Nz},
      tree(create_particle_octree(cells))
//      update_interval(),
    {
      for (auto& cell : cells) {
        cell.reserve(max_ppc);
      }
      create_particle_octree(cells);
    }

//    [[nodiscard]] bool update_this_step(const compute_t time) {
//      return update_interval.update_this_step(time);
//    }

    void add_particle(Particle&& p, const std::size_t cid) {
      cells[cid].push_back(std::move(p));
    }

    static bool update_tree_nodes(auto& node) {
      for (std::size_t i = 0; i < 8; i++) {
        if (node.is_leaf) {
          // check for particles in each cell and set the appropriate active bit
          node.active.set(i, !node.cells[i]->empty());
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
    std::size_t atomic_number;
    compute_t mass;
    compute_t charge;
    compute_t inv_mass;
    compute_t cm_ratio;
    compute_t inv_cm_ratio;
    compute_t qdt_over_2m;

//    Interval update_interval{};
    std::vector<std::vector<Particle>> cells{Nx * Ny * Nz};
    Octree<std::vector<Particle>> tree;
  }; // end struct ParticleGroup

  // struct ParticleInitializer {
  //   static void initializeFromFile(ParticleGroup& g, const std::string& filename) {
  //     std::fstream file(filename, std::ios::in);
  //
  //     if (!file.is_open()) {
  //       throw std::runtime_error("Particle initialization from file failed: " + filename);
  //     }
  //
  //     vec3<compute_t> location{};
  //     vec3<compute_t> velocity{};
  //     float weight = 0.0;
  //
  //     std::println("Opened particle file: {}", filename);
  //     std::string line;
  //     while (getline(file, line)) {
  //       std::istringstream buffer(line);
  //       buffer >> location >> velocity >> weight;// >> uid;
  //
  //       //std::cout << uid << group.mass << group.charge << std::endl;
  //
  //       // compute Lorentz factor and relativistic momentum
  //       double gamma = 1.0 / std::sqrt(1.0 - velocity.length_squared() / constants::c_sqr);
  //       auto momentum = gamma * g.mass * velocity;
  //
  //       // add particle to group
  //       g.add_particle({location, location, momentum, gamma, weight});
  //     }
  //     file.close();
  //   } // end initializeFromFile
  // }; // end struct Particle Initializer
} // end namespace tf::particles

#endif //PARTICLE_HPP
