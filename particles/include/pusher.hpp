#ifndef PUSHER_HPP
#define PUSHER_HPP

#include "particles.hpp"
#include "em_data.hpp"
#include "current_deposition.hpp"

#include "octree.hpp"

#include <cmath>

namespace tf::particles {

  struct BorisPush {
    using emdata_t = electromagnetics::EMData;
    using group_t = ParticleGroup;
    using p_tree = Octree<ParticleCell>;
    using efield = std::array<compute_t, 4>;
    using bfield = std::array<compute_t, 2>;

    static void update_particle(Particle&,
                                const efield&, const efield&, const efield&,
                                const bfield&, const bfield&, const bfield&);

    static void advance_particles(ParticleCell&, const group_t&, const emdata_t&);
    static void push_particles(const p_tree&, const group_t&, const emdata_t&);

    static void update_cells(const p_tree&, group_t&);

    static void operator()(group_t& g, const emdata_t& emdata) {
      push_particles(g.tree, g, emdata);
      update_cells(g.tree, g);
      g.update_tree();
    }
  }; // end struct BorisPush

} // end namespace tf::particles



#endif //PUSHER_HPP
