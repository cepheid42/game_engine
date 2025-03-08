#ifndef EM_ARRAY_HPP
#define EM_ARRAY_HPP

//#include <cstdint>
#include <array>

namespace tf::electromagnetics {
  template<typename T, std::size_t NX, std::size_t NY, std::size_t NZ>
  class Array3D {
  public:
    explicit Array3D() = default;
    explicit Array3D(const T fill_) : data{} { fill(fill_); }

    auto begin() { return data.begin(); }
    auto end() { return data.end(); }

    [[nodiscard]] constexpr std::size_t size() const { return data.size(); }
    [[nodiscard]] constexpr std::size_t num_bytes() const { return size() * sizeof(T); }
    [[nodiscard]] std::size_t get_scid(const std::size_t i, const std::size_t j, const std::size_t k) const { return k + nz_ * j + nz_ * ny_ * i; }

    // Specialized accessors
    T& operator()(const std::size_t i, const std::size_t j, const std::size_t k) { return (*this)[get_scid(i, j, k)]; }
    const T& operator()(const std::size_t i, const std::size_t j, const std::size_t k) const { return (*this)[get_scid(i, j, k)]; }

    T& operator[](const std::size_t i) { return data[i]; }
    const T& operator[](const std::size_t i) const { return data[i]; }

    void fill(T value) { for (auto& el : this->data) el = value; }

    // Dims
    [[nodiscard]] constexpr std::array<std::size_t, 3> dims() const { return {nx_, ny_, nz_}; }
    [[nodiscard]] constexpr std::size_t nx() const { return nx_; }
    [[nodiscard]] constexpr std::size_t ny() const { return ny_; }
    [[nodiscard]] constexpr std::size_t nz() const { return nz_; }

    // Unary Negation
    auto operator-() const;
    
    // Augmented Assignment
    auto& operator+=(const Array3D& other) {
      assert(nx_ == other.nx() && ny_ == other.ny() && nz_ == other.nz());
      for (std::size_t i = 0; i < this->data.size(); i++) { (*this)[i] += other.data[i]; }
      return *this;
    }
    
    auto& operator-=(const Array3D& other) {
      assert(nx_ == other.nx() && ny_ == other.ny() && nz_ == other.nz());
      for (std::size_t i = 0; i < this->data.size(); i++) { (*this)[i] -= other.data[i]; }
      return *this;
    }
    
    // Array-Scalar Operators
    auto& operator+=(const T s) {
      for (auto& value : this->data) { value += s; }
      return *this;
    }
    
    auto& operator-=(const T s) {
      for (auto& value : this->data) { value -= s; }
      return *this;
    }
    
    auto& operator*=(const T s) {
      for (auto& value : this->data) { value *= s; }
      return *this;
    }
    
    auto& operator/=(const T s) {
      for (auto& value : this->data) { value /= s; }
      return *this;
    }

  private:
    // Stride data
    const std::size_t nx_ = NX;
    const std::size_t ny_ = NY;
    const std::size_t nz_ = NZ;
    std::array<T, NX * NY * NZ> data;
    //
  };

}

#endif //EM_ARRAY_HPP
