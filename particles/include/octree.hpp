#ifndef OCTREE_HPP
#define OCTREE_HPP

#include "program_params.hpp"
#include "vec3.hpp"

#include <array>
#include <vector>
#include <print>
#include <bitset>

namespace tf::particles {
  template<typename T>
  T morton_encode(T x, T y, T z) {
    static constexpr std::array<T, 4> S = {16, 8, 4, 2};
    static constexpr std::array<T, 4> B = {0xFF0000FF, 0x0300F00F, 0x030C30C3, 0x09249249};

    for (std::size_t i = 0; i < 4; i++) {
      x = (x | (x << S[i])) & B[i];
      y = (y | (y << S[i])) & B[i];
      z = (z | (z << S[i])) & B[i];
    }

    return z | (y << 1) | (x << 2);
  }

  template<typename T>
  std::array<std::size_t, 3> morton_decode(T idx) {
    static constexpr std::array<T, 4> S = {2, 4, 8, 16};
    static constexpr std::array<T, 5> B = {0x030C30C3, 0x0300F00F, 0xFF0000FF, 0x000003FF};

    T x = (idx >> 2) & 0x09249249;
    T y = (idx >> 1) & 0x09249249;
    T z = idx & 0x09249249;

    for (std::size_t i = 0; i < 4; i++) {
      x = (x | (x >> S[i])) & B[i];
      y = (y | (y >> S[i])) & B[i];
      z = (z | (z >> S[i])) & B[i];
    }

    return {x, y, z};
  }

  inline auto getOctreeDepth(const std::size_t n) {
    const auto den = std::log(8.0);
    const auto num = std::log(static_cast<double>(n));
    return static_cast<std::size_t>(std::floor(num / den));
  }

  template<typename P>
  struct Octree {
    std::array<std::size_t, 3> cell_coords{};
    std::array<P*, 8> cells{}; // P = std::vector<particle>
    std::vector<Octree> children{};
    std::bitset<8> active{};
    bool is_leaf{false};
  };

  template<typename P>
  Octree<P> build_octree(std::vector<P>& cells, std::size_t begin, std::size_t end, std::size_t depth_limit) {
    Octree<P> oct{};

    if (depth_limit == 0) {
      // leaf nodes
      oct.cell_coords = morton_decode(begin);
      oct.is_leaf = true;
      for (std::size_t i = 0; i < 8; i++) {
        oct.cells[i] = &(cells[begin + i]);
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
  Octree<P> create_particle_octree(std::vector<P>& cells) {
    const auto depth = getOctreeDepth(cells.size());
    return build_octree<P>(cells, 0, cells.size(), depth - 1);
  }
}
#endif //OCTREE_HPP
