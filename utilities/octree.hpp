#ifndef OCTREE_HPP
#define OCTREE_HPP

#include "morton.hpp"
#include "vec3.hpp"

#include <vector>
#include <cstdint>
#include <cstring>
#include <cassert>
#include <cmath>
#include <bitset>
#include <numeric>


namespace tf {
  inline auto getOctreeDepth(const std::size_t n) {
    const auto den = std::log(8.0);
    const auto num = std::log(static_cast<double>(n));
    return static_cast<std::size_t>(std::floor(num / den));
  }

  inline auto nearestPowerOf2(std::size_t v) {
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v++;
    return v;
  }

  template<typename P>
  struct OctNode {
    std::array<P*, 8> cells{}; // P = std::vector<particle>
    std::vector<OctNode> children{};

    std::bitset<8> active{};
    bool is_leaf{false};
  };

  // template<typename P>
  // void visit_octree(Octree<P>& node) {
  //   for (std::size_t c = 0; c < 8; c++) {
  //     if (node.is_leaf) {
  //       if (node.cells[c]->index == 32) {
  //         node.cells[c]->value = 9999;
  //       }
  //       // const auto p = *(node.cells[c]);
  //       //std::println("{}: {}", p.index, p.value);
  //     } else {
  //       visit_octree(node.children[c]);
  //     }
  //   }
  // }

  template<typename P>
  OctNode<P> build_octree(std::vector<P>& cells, const std::size_t begin, const std::size_t end, const std::size_t depth_limit, const std::array<std::size_t, 3>& dims) {
    OctNode<P> oct{};

    if (depth_limit == 0) {
      // leaf nodes
      oct.is_leaf = true;
      for (std::size_t i = 0; i < 8; i++) {
        oct.cells[i] = &cells[begin + i];
      }
    } else {
      const auto chunk = (end - begin) / 8;
      for (std::size_t i = begin; i < end; i += chunk) {
        oct.children.push_back(build_octree<P>(cells, i, i + chunk, depth_limit - 1));
      }
    }
    return oct;
  }

  template<typename P>
  OctNode<P> create_octree(std::vector<P>& cells, const std::size_t nx, std::size_t ny, std::size_t nz) {
    const auto nearest_nx = nearestPowerOf2(nx);
    const auto nearest_ny = nearestPowerOf2(ny);
    const auto nearest_nz = nearestPowerOf2(nz);

    const auto depth = getOctreeDepth(nearest_nx * nearest_ny * nearest_nz);
    return build_octree<P>(cells, 0, cells.size(), depth - 1);
  }
} // end namespace tf
#endif //OCTREE_HPP
