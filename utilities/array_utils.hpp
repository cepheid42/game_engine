#ifndef GAME_ENGINE_ARRAY_UTILS_HPP
#define GAME_ENGINE_ARRAY_UTILS_HPP

#include "mdspan.hpp"

#include <print>

// template<typename T>
// struct Slice {
//    T mdspan;
//    std::size_t x{0};
//    std::size_t y{0};
//    std::size_t z{0};
// };
//
// constexpr auto slice1D(const auto sl) {
//    constexpr auto extent0 = sl.mdspan.extent(0);
//    constexpr auto extent1 = sl.mdspan.extent(1);
//
//    if constexpr (sl.x != 0zu) {
//       std::extents shape{ extent1 };
//       std::array strides{ extent0 * extent1 };
//       return std::mdspan{ &sl.mdspan[sl.x, 0, 0], {shape, strides} };
//    }
//    else if constexpr (sl.y != 0zu) {
//       std::extents shape{ extent0, extent2 };
//       std::array strides{ extent1 * extent2, 1zu };
//       return std::mdspan{ &sl.mdspan[0, sl.y, 0], {shape, strides} };
//    }
//    else if constexpr (sl.z != 0zu) {
//       std::extents shape{ extent0, extent1 };
//       std::array strides{ extent1 * extent2, extent2 };
//       return std::mdspan{ &sl.mdspan[0, 0, sl.z], {shape, strides} };
//    }
//    else {
//       assert(false);
//    }
// }
//
// template <typename Extents>
// requires (Extents::rank() == 3)
// constexpr auto slice(auto sl.mdspan, const auto slice_index, const auto dimension_index)
// {
//    constexpr auto extent0 = sl.mdspan.extent(0);
//    constexpr auto extent1 = sl.mdspan.extent(1);
//    constexpr auto extent2 = sl.mdspan.extent(2);
//
//    if constexpr (dimension_index == 0) {
//          std::extents shape{ extent1, extent2 };
//          std::array strides{ extent2, 1uz };
//          return std::mdspan{ &sl.mdspan[slice_index, 0, 0], {shape, strides} };
//    }
//    else if constexpr (dimension_index == 1) {
//          std::extents shape{ extent0, extent2 };
//          std::array strides{ extent1 * extent2, 1uz };
//          return std::mdspan{ &sl.mdspan[0, slice_index, 0], {shape, strides} };
//    }
//    else if constexpr (dimension_index == 2) {
//          std::extents shape{ extent0, extent1 };
//          std::array strides{ extent1 * extent2, extent2 };
//          return std::mdspan{ &sl.mdspan[0, 0, slice_index], {shape, strides} };
//    }
//    else {
//       assert(false);
//    }
// }

auto fill_mdarray(auto& A, const auto idx) {
   for (auto i = 0u; i < A.extent(0); i++) {
      for (auto j = 0u; j < A.extent(1); j++) {
         for (auto k = 0u; k < A.extent(2); k++) {
            A[i, j, k] = idx == 0 ? i : (idx == 1 ? j : k);
         }
      }
   }
}


auto diff_mdarray(auto& A, const auto& B, const auto idx) {
   for (auto i = 0u; i < A.extent(0); i++) {
      for (auto j = 0u; j < A.extent(1); j++) {
         for (auto k = 0u; k < A.extent(2); k++) {
            double val{};
            if      (idx == 0) { val = B[i + 1, j, k] - B[i, j, k]; }
            else if (idx == 1) { val = B[i, j + 1, k] - B[i, j, k]; }
            else if (idx == 2) { val = B[i, j, k + 1] - B[i, j, k]; }
            A[i, j, k] = val;
         }
      }
   }
}

auto print_mdarray(const auto& A) {
   for (auto i = 0u; i < A.extent(0); i++) {
      for (auto j = 0u; j < A.extent(1); j++) {
         for (auto k = 0u; k < A.extent(2); k++) {
            std::print("{},", A[i, j, k]);
         }
         std::println();
      }
      std::println();
   }
}

#endif //GAME_ENGINE_ARRAY_UTILS_HPP