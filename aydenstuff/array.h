//
// Created by akis on 9/4/24.
//

#ifndef TRIFORCE_ARRAY_H
#define TRIFORCE_ARRAY_H

#include <vector>

#include "tags.h"

using std::size_t;

// ===== Array =====
// =================
namespace tf::types
{
  // ----- Detail -----
  namespace detail
  {
    template <typename T, size_t N>
    class ArrayBase {
    public:
      using value_t = T;
      using vector_t = std::vector<value_t>;
      using dimension_t = tags::Dimension<N>;
      
      explicit ArrayBase() = default;
      explicit ArrayBase(size_t size) : data(size) {}
      ArrayBase(size_t size, value_t fill) : data(size, fill) {}
      ArrayBase(const ArrayBase& other) : data(other.data) {}
      
      ~ArrayBase() = default;
      
      [[nodiscard]] size_t size() const { return data.size(); }
      [[nodiscard]] size_t num_bytes() const { return size() * sizeof(value_t); }
      
      explicit operator bool() const { return data != nullptr; }
      
      auto begin() const { return data.begin(); }
      auto end() const { return data.end(); }
      
      // Basic accessors
      value_t& operator[](std::size_t i) { return data[i]; }
      const value_t& operator[](std::size_t i) const { return data[i]; }
    
    protected:
      vector_t data;
    };
  }
  
  // ----- Specialization -----


  template <typename T>
  struct Array1D : public detail::ArrayBase<T, 1> {
    using value_t = typename detail::ArrayBase<T, 1>::value_t;
    using vector_t = typename detail::ArrayBase<T, 1>::vector_t;
    using dimension_t = typename detail::ArrayBase<T, 1>::dimension_t;
    
    // Constructors & Destructor
    explicit Array1D() = default;
    explicit Array1D(std::size_t nx_) : detail::ArrayBase<T, 1>(nx_), nx(nx_) {}
    Array1D(std::size_t nx_, value_t fill) : detail::ArrayBase<T, 1>(nx_, fill), nx(nx_) {}
    Array1D(const Array1D& other) : detail::ArrayBase<T, 1>(other), nx(other.nx) {}
    
    // Array1D& operator=(const Array1D& other);
    
    ~Array1D() = default;
    
    // Indexing Function
    [[nodiscard]] static inline size_t get_scid(std::size_t i) { return i; }
    
    // Specialized Accessors
    value_t& operator()(std::size_t i) { return (*this)[i]; }
    const value_t& operator()(std::size_t i) const { return (*this)[i]; }
    
    // Unary Negation
    auto operator-() const {
      Array1D result(nx);
      for (std::size_t i = 0; i < nx; i++) { result[i] = value_t(-1) * this->data[i]; }
      return result;
    }
    
    // Augmented Assignment
    auto& operator+=(const Array1D& other) {
      assert(nx == other.nx);
      for (size_t i = 0; i < nx; i++) { this->data[i] += other.data[i]; }
      return *this;
    }
    
    auto& operator-=(const Array1D& other) {
      assert(nx == other.nx);
      for (size_t i = 0; i < nx; i++) { this->data[i] -= other.data[i]; }
      return *this;
    }
    
    // Array-Scalar Operators
    auto& operator+=(const value_t s) {
      for (auto& value : this->data) { value += s; }
      return *this;
    }
    
    auto& operator-=(const value_t s) {
      for (auto& value : this->data) { value -= s; }
      return *this;
    }
    
    auto& operator*=(const value_t s) {
      for (auto& value : this->data) { value *= s; }
      return *this;
    }
    
    auto& operator/=(const value_t s) {
      for (auto& value : this->data) { value /= s; }
      return *this;
    }
    
    // Stride data
    const size_t nx;
    //
  };// end class Array1D
  
  template <typename T>
  class Array2D : public detail::ArrayBase<T, 2> {
  public:
    using value_t = typename detail::ArrayBase<T, 2>::value_t;
    using vector_t = typename detail::ArrayBase<T, 2>::vector_t;
    using dimension_t = typename detail::ArrayBase<T, 2>::dimension_t;
    
    // Constructors & Destructor
    explicit Array2D() = default;
    Array2D(std::size_t nx_, std::size_t nz_) : detail::ArrayBase<T, 2>(nx_ * nz_), nx(nx_), nz(nz_) {}
    Array2D(std::size_t nx_, std::size_t nz_, value_t fill) : detail::ArrayBase<T, 2>(nx_ * nz_, fill), nx(nx_), nz(nz_) {}
    Array2D(const Array2D& other) : detail::ArrayBase<T, 2>(other), nx(other.nx), nz(other.nz) {}
    
    ~Array2D() = default;
    
    // Indexing Function
    [[nodiscard]] size_t get_scid(const size_t  i, const size_t  k) const { return k + nz * i; }
    
    // Specialized accessors
    value_t& operator()(const size_t  i, const size_t  k) { return (*this)[get_scid(i, k)]; }
    const value_t& operator()(const size_t  i, const size_t  k) const { return (*this)[get_scid(i, k)]; }

    value_t forward_diff_x(const size_t i, const size_t k) const { return (*this)[get_scid(i + 1, k)] - (*this)[get_scid(i, k)]; }
    value_t backward_diff_x(const size_t i, const size_t k) const { return (*this)[get_scid(i - 1, k)] - (*this)[get_scid(i, k)]; }

    value_t forward_diff_z(const size_t i, const size_t k) const { return (*this)[get_scid(i, k + 1)] - (*this)[get_scid(i, k)]; }
    value_t backward_diff_z(const size_t i, const size_t k) const { return (*this)[get_scid(i, k - 1)] - (*this)[get_scid(i, k)]; }

    // Unary Negation
    auto operator-() const {
      Array2D result(nx);
      for (std::size_t i = 0; i < nx; i++) { result[i] = value_t(-1) * this->data[i]; }
      return result;
    }
    
    // Augmented Assignment
    auto& operator+=(const Array2D& other) {
      assert(nx == other.nx);
      for (size_t i = 0; i < nx; i++) { this->data[i] += other.data[i]; }
      return *this;
    }
    
    auto& operator-=(const Array2D& other) {
      assert(nx == other.nx);
      for (size_t i = 0; i < nx; i++) { this->data[i] -= other.data[i]; }
      return *this;
    }
    
    // Array-Scalar Operators
    auto& operator+=(const value_t s) {
      for (auto& value : this->data) { value += s; }
      return *this;
    }
    
    auto& operator-=(const value_t s) {
      for (auto& value : this->data) { value -= s; }
      return *this;
    }
    
    auto& operator*=(const value_t s) {
      for (auto& value : this->data) { value *= s; }
      return *this;
    }
    
    auto& operator/=(const value_t s) {
      for (auto& value : this->data) { value /= s; }
      return *this;
    }
    
    // Stride data
    const std::size_t nx, nz;
    //
  };// end class Array2D
  
  template <typename T>
  class Array3D : public detail::ArrayBase<T, 3> {
  public:
    using value_t = typename detail::ArrayBase<T, 3>::value_t;
    using vector_t = typename detail::ArrayBase<T, 3>::vector_t;
    using dimension_t = typename detail::ArrayBase<T, 3>::dimension_t;
    
    Array3D(size_t nx_, size_t ny_, size_t nz_) : detail::ArrayBase<T, 3>(nx_ * ny_ * nz_), nx(nx_), ny(ny_), nz(nz_) {}
    
    inline size_t get_scid(std::size_t i, std::size_t j, std::size_t k) { return k + nz * j + nz * ny * i; }
    
    // Specialized accessors
    value_t& operator()(std::size_t i, std::size_t j, std::size_t k) { return (*this)[get_scid(i, j, k)]; }
    const value_t& operator()(std::size_t i, std::size_t j, std::size_t k) const { return (*this)[get_scid(i, j, k)]; }
    
    // Stride data
    const size_t nx, ny, nz;
    //
  };// end class Array3D
  //
} // end namespace tf::types

#endif //TRIFORCE_ARRAY_H
