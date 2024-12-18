//
// Created by cepheid on 10/21/24.
//

#ifndef CURL_OPERATORS_H
#define CURL_OPERATORS_H

namespace tf::electromagnetics
{
  using tf::electromagnetics::traits::Derivative;

  //=================== Array Differencing Functions ========================
  //=========================================================================
  // NoOp/Default Diff function
  template<Derivative, bool>
  struct Diff {
    static constexpr auto apply(const auto&, size_t) { return 0.0; }
    static constexpr auto apply(const auto&, size_t, size_t) { return 0.0; }
    static constexpr auto apply(const auto&, size_t, size_t, size_t) { return 0.0; }
  };

  // ================= Forward Differences =================
  template<>
  struct Diff<Derivative::DX, true> {
    static auto apply(const auto& f, size_t i)                     { return f(i + 1) - f(i); }
    static auto apply(const auto& f, size_t i, size_t j)           { return f(i + 1, j) - f(i, j); }
    static auto apply(const auto& f, size_t i, size_t j, size_t k) { return f(i + 1, j, k) - f(i, j, k); }
  };

  template<>
  struct Diff<Derivative::DY, true> {
    static auto apply(const auto& f, size_t i, size_t j)           { return f(i, j + 1) - f(i, j); }
    static auto apply(const auto& f, size_t i, size_t j, size_t k) { return f(i, j + 1, k) - f(i, j, k); }
  };

  template<>
  struct Diff<Derivative::DZ, true> {
    static auto apply(const auto& f, size_t i, size_t j, size_t k) { return f(i, j, k + 1) - f(i, j, k); }
  };

  // ================= Backward Differences =================
  template<>
  struct Diff<Derivative::DX, false> {
    static auto apply(const auto& f, size_t i)                     { return f(i) - f(i - 1); }
    static auto apply(const auto& f, size_t i, size_t j)           { return f(i, j) - f(i - 1, j); }
    static auto apply(const auto& f, size_t i, size_t j, size_t k) { return f(i, j, k) - f(i - 1, j, k); }
  };

  template<>
  struct Diff<Derivative::DY, false> {
    static auto apply(const auto& f, size_t i, size_t j)           { return f(i, j) - f(i, j - 1); }
    static auto apply(const auto& f, size_t i, size_t j, size_t k) { return f(i, j, k) - f(i, j - 1, k); }
  };
  template<>
  struct Diff<Derivative::DZ, false> {
    static auto apply(const auto& f, size_t i, size_t j, size_t k) { return f(i, j, k) - f(i, j, k - 1); }
  };
}
#endif //CURL_OPERATORS_H
