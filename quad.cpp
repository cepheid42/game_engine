// #include "/home/cepheid/TriForce/game_engine/utilities/morton.hpp"

#include <print>
#include <cstdint>
#include <array>
#include <vector>
#include <bitset>
#include <bit>
#include <cmath>
#include <algorithm>
#include <utility>
#include <cassert>


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

struct Point {
  std::size_t x, y;
};

struct CellData {
  Point idx{};
  std::uint32_t code{};
  std::vector<int> chunks{};
};

struct Box {
  Point min;
  Point max;
};

inline constexpr auto null = static_cast<std::uint32_t>(-1);

struct Node {
  CellData* cell{};
  std::array<std::uint32_t, 4> children{};
  std::bitset<4> active;
  bool is_leaf{};
};

struct Quadtree {
  std::uint32_t root{};
  std::vector<Node> nodes{};
};

template<typename Iterator>
std::uint32_t build_impl(Quadtree& tree, const Box& box, Iterator begin, Iterator end, const std::size_t depth) {
  if (begin == end) { assert(false); }

  const auto result = tree.nodes.size();
  tree.nodes.emplace_back();

  if (begin + 1 == end or depth == 0) {
    tree.nodes[result].is_leaf = true;
    tree.nodes[result].cell = &(*begin);
    return result;
  }

  auto middle = [](const Point& a, const Point& b) -> Point {
    return {(a.x + b.x) / 2, (a.y + b.y) / 2};
  };

  const auto center = middle(box.min, box.max);

  auto bottom = [&](const CellData& p) { return p.idx.y < center.y; };
  auto left = [&](const CellData& p) { return p.idx.x < center.x; };

  Iterator split_y = std::partition(begin, end, bottom);
  Iterator split_x_lower = std::partition(begin, split_y, left);
  Iterator split_x_upper = std::partition(split_y, end, left);

  tree.nodes[result].children[0] = build_impl(tree, {box.min, center}, begin, split_x_lower, depth - 1);
  tree.nodes[result].children[1] = build_impl(tree, {{center.x, box.min.y}, {box.max.x, center.y}},split_x_lower, split_y, depth - 1);
  tree.nodes[result].children[2] = build_impl(tree, {{box.min.x, center.y}, {center.x, box.max.y}}, split_y, split_x_upper, depth - 1);
  tree.nodes[result].children[3] = build_impl(tree, {center, box.max}, split_x_upper, end, depth - 1);

  return result;
}

Quadtree build_quadtree(std::vector<CellData>& cells) {
  constexpr auto tree_depth = 10;
  Quadtree result{};
  result.root = build_impl(result, Box{cells.begin()->idx, cells.end()->idx}, cells.begin(), cells.end(), tree_depth);
  return result;
}

void visit(const Node& node, const Quadtree& tree) {
  if (node.is_leaf) {
    std::println("{}, {}", node.cell->idx.x, node.cell->idx.y);
  } else {
    for (const auto& c : node.children) {
      visit(tree.nodes[c], tree);
    }
  }
}


int main() {
  constexpr std::size_t nx = 8, ny = 8, nz = 8;

  std::vector<CellData> cells;

  for (std::size_t i = 0; i < nx; i++) {
    for (std::size_t j = 0; j < ny; j++) {
      const auto idx = j + ny * i;
      cells.push_back({.idx = {i, j}});
    }
  }
  auto tree = build_quadtree(cells);

  visit(tree.nodes[tree.root], tree);


  return 0;
}