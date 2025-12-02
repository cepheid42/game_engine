#ifndef GAME_ENGINE_ARRAY_UTILS_HPP
#define GAME_ENGINE_ARRAY_UTILS_HPP

#include "mdspan.hpp"

#include <print>

template<typename mdspan_t>
constexpr auto strided_slice(const mdspan_t& mdspan) {}

template<typename mdspan_t>
requires (mdspan_t::rank() == 3)
constexpr auto slice3D(const mdspan_t& mdspan, const std::array<std::size_t, 3>& idxs, const std::dextents<std::size_t, 3>& ext)
-> std::mdspan<typename mdspan_t::value_type, std::dextents<std::size_t, 3>, std::layout_stride>
{
   const auto& [i, j, k] = idxs;
   const auto& nx = mdspan.extent(0);
   const auto& ny = mdspan.extent(1);
   const auto& nz = mdspan.extent(2);
   return {&mdspan[i, j, k], {ext, std::array{ny * nz, nz, 1zu}}};
}

template<typename mdspan_t>
requires (mdspan_t::rank() == 3)
constexpr auto slice2D(const mdspan_t& mdspan, const auto idx, const auto dim)
-> std::mdspan<typename mdspan_t::value_type, std::dextents<std::size_t, 2>, std::layout_stride>
{
   const auto nx = mdspan.extent(0);
   const auto ny = mdspan.extent(1);
   const auto nz = mdspan.extent(2);
   if (dim == 0) {
      return {&mdspan[idx, 0, 0], {std::extents{ny, nz}, std::array{nz, 1zu}}};
   }
   if (dim == 1) {
      return {&mdspan[0, idx, 0], {std::extents{nx, nz}, std::array{ny * nz, 1zu}}};
   }
   if (dim == 2) {
      return {&mdspan[0, 0, idx], {std::extents{nx, ny}, std::array{ny * nz, nz}}};
   }
   assert(false);
}

template<typename mdspan_t>
requires (mdspan_t::rank() == 3)
constexpr auto slice1D(const mdspan_t& mdspan, const auto idx1, const auto idx2, const auto dim)
-> std::mdspan<typename mdspan_t::value_type, std::dextents<std::size_t, 1>, std::layout_stride>
{
   const auto nx = mdspan.extent(0);
   const auto ny = mdspan.extent(1);
   const auto nz = mdspan.extent(2);
   if (dim == 0) {
      return mdspan_t{&mdspan[0, idx1, idx2], {std::extents{nx}, std::array{ny * nz}}};
   }
   if (dim == 1) {
      return mdspan_t{&mdspan[idx1, 0, idx2], {std::extents{ny}, std::array{nz}}};
   }
   if (dim == 2) {
      return mdspan_t{&mdspan[idx1, idx2, 0], {std::extents{nz}, std::array{1zu}}};
   }
   assert(false);
}

auto stepfill_mdarray(auto& A, const auto idx) {
   for (auto i = 0u; i < A.extent(0); i++) {
      for (auto j = 0u; j < A.extent(1); j++) {
         for (auto k = 0u; k < A.extent(2); k++) {
            A[i, j, k] = idx == 0 ? i : (idx == 1 ? j : k);
         }
      }
   }
}

template<typename mdspan_t>
auto fill_mdarray(mdspan_t& A, const auto val) requires (mdspan_t::rank() == 3) {
   for (auto i = 0u; i < A.extent(0); i++) {
      for (auto j = 0u; j < A.extent(1); j++) {
         for (auto k = 0u; k < A.extent(2); k++) {
            A[i, j, k] = val;
         }
      }
   }
}

template<typename mdspan_t>
auto fill_mdarray(mdspan_t& A, const auto val) requires (mdspan_t::rank() == 2) {
   for (auto i = 0u; i < A.extent(0); i++) {
      for (auto j = 0u; j < A.extent(1); j++) {
            A[i, j] = val;
      }
   }
}

template<typename mdspan_t>
auto print_mdarray(const mdspan_t& A) requires (mdspan_t::rank() == 3) {
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

template<typename mdspan_t>
auto print_mdarray(const mdspan_t& A) requires (mdspan_t::rank() == 2) {
   for (auto i = 0u; i < A.extent(0); i++) {
      for (auto j = 0u; j < A.extent(1); j++) {
         std::print("{},", A[i, j]);
      }
      std::println();
   }
}
#endif //GAME_ENGINE_ARRAY_UTILS_HPP