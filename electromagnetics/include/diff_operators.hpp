#ifndef EM_CURL_HPP
#define EM_CURL_HPP

namespace tf::electromagnetics
{
  enum class Derivative { DX, DY, DZ };
  
  //=================== Array Differencing Functions ========================
  //=========================================================================
  template<Derivative, bool>
  struct Diff;

  // ================= Forward Differences =================
  template<Derivative D>
  struct Diff<D, true> {
    static constexpr Derivative type = D;
    static constexpr bool Forward = true;

    static auto apply(const auto& f, std::size_t i, std::size_t j, std::size_t k) {
      if constexpr (D == Derivative::DX) {
        return f(i + 1, j, k) - f(i, j, k);
      }
      else if constexpr (D == Derivative::DY) {
        return f(i, j + 1, k) - f(i, j, k);
      }
      else {
        return f(i, j, k + 1) - f(i, j, k);
      }
    }
  };

  // ================= Backward Differences =================
  template<Derivative D>
  struct Diff<D, false> {
    static constexpr Derivative type = D;
    static constexpr bool Forward = false;

    static auto apply(const auto& f, std::size_t i, std::size_t j, std::size_t k) {
      if constexpr (D == Derivative::DX) {
        return f(i, j, k) - f(i - 1, j, k);
      }
      else if constexpr (D == Derivative::DY) {
        return f(i, j, k) - f(i, j - 1, k);
      }
      else {
        return f(i, j, k) - f(i, j, k - 1);
      }
    }
  };

  using forward_dx = Diff<Derivative::DX, true>;
  using forward_dy = Diff<Derivative::DY, true>;
  using forward_dz = Diff<Derivative::DZ, true>;
  using backward_dx = Diff<Derivative::DX, false>;
  using backward_dy = Diff<Derivative::DY, false>;
  using backward_dz = Diff<Derivative::DZ, false>;
} // end namespace tf::electromagnetics

#endif //EM_CURL_HPP
