#ifndef COLLISIONS_HPP
#define COLLISIONS_HPP

#include "particles.hpp"

#include <vector>

namespace tf::collisions
{
struct Collisions {
   std::vector<std::size_t> offsets{};

   void find_idx_offsets(const auto& g) {
      offsets.emplace_back(0); // first particle
      auto old_code = morton_encode(particles::getCellIndices(g.particles[0].location)); // first particles code

      // loop over remaining particles
      for (std::size_t i = 1; i < g.particles.size(); ++i) {
         auto new_code = morton_encode(particles::getCellIndices(g.particles[i].location));
         if (old_code != new_code) {
            offsets.emplace_back(i);
            old_code = new_code;
         }
      }
      // last offset is 1 past the end
      offsets.emplace_back(g.particles.size());
   }


};
} // end namespace tf::collisions

#endif //COLLISIONS_HPP
