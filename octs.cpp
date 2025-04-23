#include "/home/cepheid/TriForce/game_engine/utilities/morton.hpp"


#include <array>
#include <vector>
#include <print>
#include <bitset>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <utility>

// inline constexpr std::array<std::uint32_t, 8> CHILD_CODES = {
//   0b000, // TLF
//   0b001, // TRF
//   0b010, // TLB
//   0b011, // TRB
//   0b100, // BLF
//   0b101, // BRF
//   0b110, // BLB
//   0b111  // BRF
// };

// template<typename T>
// T morton_encode(T x, T y, T z) {
//   static constexpr std::array<T, 4> S = {16, 8, 4, 2};
//   static constexpr std::array<T, 4> B = {0xFF0000FF, 0x0300F00F, 0x030C30C3, 0x09249249};
//
//   for (std::size_t i = 0; i < 4; i++) {
//     x = (x | (x << S[i])) & B[i];
//     y = (y | (y << S[i])) & B[i];
//     z = (z | (z << S[i])) & B[i];
//   }
//
//   return z | (y << 1) | (x << 2);
// }
//
// template<typename T>
// std::array<std::size_t, 3> morton_decode(T idx) {
//   static constexpr std::array<T, 4> S = {2, 4, 8, 16};
//   static constexpr std::array<T, 5> B = {0x030C30C3, 0x0300F00F, 0xFF0000FF, 0x000003FF};
//
//   T x = (idx >> 2) & 0x09249249;
//   T y = (idx >> 1) & 0x09249249;
//   T z = idx & 0x09249249;
//
//   for (std::size_t i = 0; i < 4; i++) {
//     x = (x | (x >> S[i])) & B[i];
//     y = (y | (y >> S[i])) & B[i];
//     z = (z | (z >> S[i])) & B[i];
//   }
//
//   return {x, y, z};
// }

inline auto getOctreeDepth(const std::size_t n) {
  const auto den = std::log(8.0);
  const auto num = std::log(static_cast<double>(n));
  return static_cast<std::size_t>(std::floor(num / den));
}

struct CellData {
  std::array<std::size_t, 3> cidxs{};
  std::vector<int> chunks{};
  std::size_t code{};
};

template<typename P>
struct Node {
  P* cell;
  std::vector<Node> children{};
  // std::bitset<8> active{};
  std::size_t code{};
  bool is_leaf{false};
};

template<typename P>
struct Octree {
  Node<P> root;
  std::size_t depth{};
};


template<typename P>
Node<P> build_octree(auto begin, auto end, std::size_t depth_limit) {
  Node<P> oct{};
  if (depth_limit == 0 or begin + 1 == end) {
    // leaf nodes
    oct.is_leaf = true;
    oct.cell = &(*begin);
    return oct;
  }

  auto start = (*begin).code;

  for (auto it = begin; it != it + 8; ++it) {
    while (start < (*it).code) { ++start; }

    std::println("{}, {}", (*it).code, start);
  }

  // const auto chunk = (end - begin) / 8;
  // for (std::size_t i = begin; i < end; i += chunk) {
  //   // std::println("{} + {} = {}", i, chunk, i + chunk);
  //   oct.children.push_back(build_octree<P>(cells, i, i + chunk, depth_limit - 1));
  // }

  return oct;
}

template<typename P>
Octree<P> create_particle_octree(std::vector<P>& cells) {
  const auto depth = getOctreeDepth(cells.size());
  Octree<P> tree;
  tree.root = build_octree<P>(cells.begin(), cells.end(), depth);
  tree.depth = depth;
  return tree;
}


// template<typename P>
// void visit_octree(Node<P>& node, const auto& Func) {
//   if (node.is_leaf) {
//     for (std::size_t i = 0; i < 8; i++) {
//       if (node.cells[i] == nullptr) { continue; }
//       std::println("Cell {}, ({}, {}, {})", node.cells[i]->code, node.cells[i]->cidxs[0], node.cells[i]->cidxs[1], node.cells[i]->cidxs[2]);
//       Func(node.cells[i]);
//     }
//   } else {
//     for (std::size_t i = 0; i < 8; i++) {
//       visit_octree(node.children[i], Func);
//     }
//   }
// }
//
// struct PrintFunctor {
//   static void operator()(const CellData* cell) {
//     for (auto& chunk : cell->chunks) {
//       if (chunk != 0) {
//         std::println("Cell {}, ({}, {}, {}) = {}", cell->code, cell->cidxs[0], cell->cidxs[1], cell->cidxs[2], chunk);
//       }
//     }
//   }
// };



int main() {
  constexpr std::size_t nx = 8;
  constexpr std::size_t ny = 1;
  constexpr std::size_t nz = 7;

  std::vector<CellData> cells;

  for (std::size_t i = 0; i < nx; i++) {
    for (std::size_t j = 0; j < ny; j++) {
      for (std::size_t k = 0; k < nz; k++) {
        const auto code = tf::morton3D_64_encode(i, j, k);
        cells.push_back({.cidxs = {i, j, k}, .chunks = {0}, .code = code});
      }
    }
  }

  std::ranges::sort(cells, {}, &CellData::code);

  // auto code000 = tf::morton3D_64_encode(0, 0, 0);
  // auto code001 = tf::morton3D_64_encode(0, 0, 1);
  // // auto code010 = tf::morton3D_64_encode(0, 1, 0);
  // auto code100 = tf::morton3D_64_encode(1, 0, 0);

  // cells[code000].chunks = {111111};
  // cells[code001].chunks = {111112};
  // // cells[code010].chunks = {111121};
  // cells[code100].chunks = {111211};

  auto tree = create_particle_octree(cells);
  //
  // visit_octree(tree.root, PrintFunctor{});
  //
  // for (std::size_t i = 0; i < nx; i++) {
  //   for (std::size_t j = 0; j < ny; j++) {
  //     for (std::size_t k = 0; k < nz; k++) {
  //       // const auto idx = tf::morton3D_64_encode(i, j, k);
  //       const auto idx = k + nz * j + ny * nz * i;
  //       // std::print("{}, ", idx);
  //       if (cells[idx].chunks[0] == 0) {
  //         std::print("{:6d}, ", cells[idx].code);
  //       } else {
  //         std::print("{:6d}, ", cells[idx].chunks[0]);
  //       }
  //     }
  //     std::println();
  //   }
  //   std::println();
  // }

  return 0;
}