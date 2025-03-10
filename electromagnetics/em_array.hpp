#ifndef EM_ARRAY_HPP
#define EM_ARRAY_HPP

#include <vector>

#include "em_vec3.hpp"

namespace tf::electromagnetics {
  template<typename T>
  class Array3D {
    using value_t = T;
  public:
    explicit Array3D() = default;
    explicit Array3D(const std::size_t nx, const std::size_t ny=1, const std::size_t nz=1)
    : nx_(nx), ny_(ny), nz_(nz), data_{nx * ny * nz} {
      std::cout << "Array3D ctor" << std::endl;
      std::cout << "nx = " << nx << ", ny = " << ny << ", nz = " << nz << std::endl;
    }

    auto begin() { return data_.begin(); }
    auto end() { return data_.end(); }

    [[nodiscard]] constexpr std::size_t size() const { return data_.size(); }
    [[nodiscard]] constexpr std::size_t num_bytes() const { return size() * sizeof(T); }
    [[nodiscard]] std::size_t get_scid(const std::size_t i, const std::size_t j, const std::size_t k) const { return k + nz_ * j + nz_ * ny_ * i; }

    // Specialized accessors
    T& operator()(const std::size_t i, const std::size_t j, const std::size_t k) { return (*this)[get_scid(i, j, k)]; }
    const T& operator()(const std::size_t i, const std::size_t j, const std::size_t k) const { return (*this)[get_scid(i, j, k)]; }

    T& operator[](const std::size_t i) { return data_[i]; }
    const T& operator[](const std::size_t i) const { return data_[i]; }

    void fill(T value) { for (auto& el : data_) el = value; }

    // Dims
    [[nodiscard]] tf::types::vec3<std::size_t> dims() const { return {nx_, ny_, nz_}; }
    [[nodiscard]] std::size_t nx() const { return nx_; }
    [[nodiscard]] std::size_t ny() const { return ny_; }
    [[nodiscard]] std::size_t nz() const { return nz_; }

    // Unary Negation
    auto operator-() const {
      Array3D<T> out(*this);
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
    std::size_t nx_;
    std::size_t ny_;
    std::size_t nz_;
    std::vector<T> data_;
  }; // end class Array3D

  template<>
  class Array3D<void> {
    using value_t = void;
  public:
    explicit Array3D() = default;
    explicit Array3D(const std::size_t, const std::size_t, const std::size_t, const auto) {}

    // Specialized accessors
    static constexpr auto operator()(const std::size_t, const std::size_t, const std::size_t) { return 0; }
    static constexpr auto operator[](const std::size_t) { return 0; }

    static constexpr void fill(auto) {}

//    // Dims
//    [[nodiscard]] tf::types::vec3<std::size_t> dims() const { return {0zu, 0zu, 0zu}; }
//    [[nodiscard]] std::size_t nx() const { return 0zu; }
//    [[nodiscard]] std::size_t ny() const { return 0zu; }
//    [[nodiscard]] std::size_t nz() const { return 0zu; }
  };
}

#endif //EM_ARRAY_HPP
