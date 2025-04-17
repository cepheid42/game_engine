#ifndef OCTREE_HPP
#define OCTREE_HPP

#include "morton.hpp"
#include "vec3.hpp"

#include <vector>
#include <cstdint>
#include <cstring>
#include <cassert>
#include <cmath>
#include <numeric>


namespace tf {
  template<typename T>
  struct CellData {
    std::uint64_t code{};
    std::vector<T> data;

    CellData() = default;
    CellData(const std::uint64_t cid, std::vector<T> d) : code(cid), data(d) {}

    bool operator<(const CellData& other) const {
      return code < other.code;
    }

    bool operator>(const CellData& other) const {
      return code > other.code;
    }
  };

  // template<typename T>
  struct Node {
    std::array<std::uint32_t, 8> children;
  };


  // template<typename T>
  struct Octree {
    std::uint32_t root;
    std::vector<Node> nodes;
  }; // end struct Octree


  std::uint32_t build_octree(Octree& tree, std::uint32_t begin, std::uint32_t end, std::size_t depth_limit) {
    if (begin == end) {
      return static_cast<std::uint32_t>(-1);
    }

    const auto result = tree.nodes.size();
    tree.nodes.emplace_back();

    if (begin + 1 == end) {
      return result;
    }

    const auto chunk = (end - begin) / 8;
    for (std::size_t i = 0; i < 8; ++i) {
      tree.nodes[result].children[i] = build_octree(tree, begin, begin + chunk, depth_limit - 1);
    }

  }

  template<typename Iterator>
  Octree create_octree(Octree& tree, Iterator begin, Iterator end) {

  }

} // end namespace tf
#endif //OCTREE_HPP
