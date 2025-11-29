#ifndef GAME_ENGINE_ARRAY_UTILS_HPP
#define GAME_ENGINE_ARRAY_UTILS_HPP

#include <print>

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