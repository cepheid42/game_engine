#ifndef PUSHER_HPP
#define PUSHER_HPP

#include "particles.hpp"
#include "em_data.hpp"

#include <cmath>

namespace tf::particles {

  struct BorisPush {
    using emdata_t = electromagnetics::EMData;
    using group_t = ParticleGroup;
    using cell_data = ParticleGroup::CellData;
    using efield = std::array<compute_t, 4>;
    using bfield = std::array<compute_t, 2>;
    using buffer_t = std::vector<std::tuple<Particle, std::size_t>>;

    buffer_t buffer{};

    static void update_velocity(Particle&,
                                const efield&, const efield&, const efield&,
                                const bfield&, const bfield&, const bfield&);
    static bool update_position(Particle&);
    static void advance_cell(cell_data&, const emdata_t&, compute_t);
    static void advance(group_t&, const emdata_t&);
    void update_cells(group_t&);

    void operator()(group_t& g, const emdata_t& emdata) {
      advance(g, emdata);
      update_cells(g);
    }
  }; // end struct BorisPush
} // end namespace tf::particles


#endif //PUSHER_HPP
