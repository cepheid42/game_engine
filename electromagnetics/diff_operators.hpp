#ifndef EM_CURL_HPP
#define EM_CURL_HPP

#include <cstddef>

namespace tf::electromagnetics {
enum class Derivative { DX, DY, DZ, NoOp };

//=================== Array Differencing Functions ========================
//=========================================================================
template<std::size_t N, Derivative, bool>
requires (N == 2)
struct Diff {
   static constexpr auto Type = Derivative::NoOp;
   static constexpr std::size_t Forward = 0;
   static auto apply(const auto&, const std::size_t, const std::size_t, const std::size_t) { return 0.0; }
};

// ================= Forward Differences =================
template<std::size_t N, Derivative D>
requires (N != 2)
struct Diff<N, D, true> {
   static constexpr Derivative  type    = D;
   static constexpr std::size_t Forward = 1;
   static auto apply(const auto& f, const std::size_t i, const std::size_t j, const std::size_t k) {
      if constexpr (D == Derivative::DX) {
         return f[i + 1, j, k] - f[i, j, k];
      }
      else if constexpr (D == Derivative::DY) {
         return f[i, j + 1, k] - f[i, j, k];
      }
      else {
         return f[i, j, k + 1] - f[i, j, k];
      }
   }
};

// ================= Backward Differences =================
template<std::size_t N, Derivative D>
requires (N != 2)
struct Diff<N, D, false> {
   static constexpr Derivative  type    = D;
   static constexpr std::size_t Forward = 0;
   static auto apply(const auto& f, const std::size_t i, const std::size_t j, const std::size_t k) {
      if constexpr (D == Derivative::DX) {
         return f[i, j, k] - f[i - 1, j, k];
      }
      else if constexpr (D == Derivative::DY) {
         return f[i, j, k] - f[i, j - 1, k];
      }
      else {
         return f[i, j, k] - f[i, j, k - 1];
      }
   }
};



// using noop        = Diff<Derivative::NoOp, true>;
template<std::size_t N> using forward_dx  = Diff<N, Derivative::DX, true>;
template<std::size_t N> using forward_dy  = Diff<N, Derivative::DY, true>;
template<std::size_t N> using forward_dz  = Diff<N, Derivative::DZ, true>;
template<std::size_t N> using backward_dx = Diff<N, Derivative::DX, false>;
template<std::size_t N> using backward_dy = Diff<N, Derivative::DY, false>;
template<std::size_t N> using backward_dz = Diff<N, Derivative::DZ, false>;
} // end namespace tf::electromagnetics

#endif //EM_CURL_HPP
