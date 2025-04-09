#ifndef OCTREE_HPP
#define OCTREE_HPP

#include "morton.hpp"
#include "vec3.hpp"

#include <vector>
#include <cstdint>
#include <cstring>
#include <cassert>
#include <cmath>


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

  template<typename T>
  struct Node {

  };


  template<typename T>
  struct Octree {

  }; // end struct Octree
} // end namespace tf
#endif //OCTREE_HPP
