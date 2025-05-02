#ifndef PUSHER_HPP
#define PUSHER_HPP

#include "em_params.hpp"
#include "particles.hpp"
#include "em_data.hpp"

#include <cmath>

namespace tf::particles {

  struct BorisPush {
    using emdata_t = electromagnetics::EMData;
    using group_t = ParticleGroup;
    using cell_data = ParticleGroup::CellData;
    using buffer_t = std::vector<std::tuple<Particle, std::size_t>>;

    static constexpr std::size_t BC_DEPTH = PMLDepth - 5;
    buffer_t buffer{};

    static void update_velocity(Particle&,
                                const auto&, const auto&, const auto&,
                                const auto&, const auto&, const auto&, auto);
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
