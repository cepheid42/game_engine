#ifndef PUSHER_HPP
#define PUSHER_HPP

#include "particles.hpp"
#include "em_data.hpp"

#include <cmath>

namespace tf::particles {

  struct BorisPush {
    using emdata_t = electromagnetics::EMData;
    using group_t = ParticleGroup;
    using p_vec = std::vector<Particle>;
    using p_tree = Octree<p_vec>;
    using efield_comp = std::array<compute_t, 4>;
    using bfield_comp = std::array<compute_t, 2>;

    static void update_particle(Particle&, const ParticleGroup&,
                                const efield_comp&, const efield_comp&, const efield_comp&,
                                const bfield_comp&, const bfield_comp&, const bfield_comp&);

    static void update_cell(p_vec&, const std::array<std::size_t, 3>&, const group_t&, const emdata_t&);

    static void visit(const p_tree&, const group_t&, const emdata_t&);

    static void operator()(const group_t& g, const emdata_t& emdata) {
       visit(g.tree, g, emdata);
    }
  }; // end struct MomentumIntegrator
} // end namespace tf::particles



#endif //PUSHER_HPP
