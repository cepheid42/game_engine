#ifndef EM_ARRAY_HPP
#define EM_ARRAY_HPP

#include "vec3.hpp"

#include <vector>

namespace tf {
  template<typename T>
  class Array3D {
    using value_t = T;
  public:
    explicit Array3D() = default;
    explicit Array3D(const std::size_t nx, const std::size_t ny, const std::size_t nz)
    : nx_(nx), ny_(ny), nz_(nz), data_(nx * ny * nz)
    {}

    auto begin() { return data_.begin(); }
    auto end() { return data_.end(); }

    auto* data() { return data_.data(); }
    [[nodiscard]] auto capacity() const { return data_.capacity(); }
    [[nodiscard]] constexpr std::size_t size() const { return data_.size(); }
    [[nodiscard]] constexpr std::size_t num_bytes() const { return size() * sizeof(T); }
    [[nodiscard]] std::size_t get_scid(const std::size_t i, const std::size_t j, const std::size_t k) const { return k + nz_ * j + nz_ * ny_ * i; }

    // Specialized accessors
    T& operator()(const std::size_t i, const std::size_t j, const std::size_t k) { return data_[get_scid(i, j, k)]; }
    const T& operator()(const std::size_t i, const std::size_t j, const std::size_t k) const { return data_[get_scid(i, j, k)]; }

    T& operator[](const std::size_t i) { return data_[i]; }
    const T& operator[](const std::size_t i) const { return data_[i]; }

    void fill(T value) { for (auto& el : data_) el = value; }

    // Dims
    [[nodiscard]] vec3<std::size_t> dims() const { return {nx_, ny_, nz_}; }
    [[nodiscard]] std::size_t nx() const { return nx_; }
    [[nodiscard]] std::size_t ny() const { return ny_; }
    [[nodiscard]] std::size_t nz() const { return nz_; }

    // Unary Negation
    auto operator-() const {
      Array3D out(*this);
      for (std::size_t i = 0; i < data_.size(); i++) { out[i] = -data_[i]; }
      return out;
    }
    
    // Augmented Assignment
    auto& operator+=(const Array3D& other) {
      assert(nx_ == other.nx() && ny_ == other.ny() && nz_ == other.nz());
      for (std::size_t i = 0; i < data_.size(); i++) { (*this)[i] += other.data_[i]; }
      return *this;
    }
    
    auto& operator-=(const Array3D& other) {
      assert(nx_ == other.nx() && ny_ == other.ny() && nz_ == other.nz());
      for (std::size_t i = 0; i < data_.size(); i++) { (*this)[i] -= other.data_[i]; }
      return *this;
    }
    
    // Array-Scalar Operators
    auto& operator+=(const T s) {
      for (auto& value : data_) { value += s; }
      return *this;
    }
    
    auto& operator-=(const T s) {
      for (auto& value : data_) { value -= s; }
      return *this;
    }
    
    auto& operator*=(const T s) {
      for (auto& value : data_) { value *= s; }
      return *this;
    }
    
    auto& operator/=(const T s) {
      for (auto& value : data_) { value /= s; }
      return *this;
    }

  private:
    // Stride data_
    std::size_t nx_{};
    std::size_t ny_{};
    std::size_t nz_{};
    std::vector<T> data_{};
  }; // end class Array3D

  template<>
  class Array3D<void> {
    using value_t = void;
  public:
    explicit Array3D() = default;
    explicit Array3D(const std::size_t, const std::size_t, const std::size_t, const auto) {}

    // Specialized accessors
    static constexpr auto operator()(const std::size_t, const std::size_t, const std::size_t) { return 0.0f; }
    static constexpr auto operator[](const std::size_t) { return 0.0f; }

    static constexpr void fill(auto) {}
  };
} // end namespace tf

#endif //EM_ARRAY_HPP
