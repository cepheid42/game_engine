//
// Created by akis on 9/4/24.
//

#ifndef TRIFORCE_ARRAY_H
#define TRIFORCE_ARRAY_H

#include <cassert>
#include <vector>

#include "tags.h"
#include "vector.h"

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
      ArrayBase(size_t size) : data(size) {}
      ArrayBase(size_t size, value_t fill) : data(size, fill) {}
      ArrayBase(const ArrayBase& other) : data(other.data) {}
      
      ~ArrayBase() = default;
      
      [[nodiscard]] size_t size() const { return data.size(); }
      [[nodiscard]] size_t num_bytes() const { return size() * sizeof(value_t); }
      
      explicit operator bool() const { return data != nullptr; }
      
      auto begin() const { return data.begin(); }
      auto end() const { return data.end(); }
      
      // Basic accessors
      value_t& operator[](size_t i) { return data[i]; }
      const value_t& operator[](size_t i) const { return data[i]; }
    
    protected:
      vector_t data;
    };
  }
  
  // ----- Array1D -----
  template <typename T>
  struct Array1D : public detail::ArrayBase<T, 1> {
    using value_t = typename detail::ArrayBase<T, 1>::value_t;
    using vector_t = typename detail::ArrayBase<T, 1>::vector_t;
    using dimension_t = typename detail::ArrayBase<T, 1>::dimension_t;
    
    // Constructors & Destructor
    explicit Array1D() = default;
    Array1D(size_t nx_) : detail::ArrayBase<T, 1>(nx_), nx_(nx_) {}
    Array1D(size_t nx_, value_t fill) : detail::ArrayBase<T, 1>(nx_, fill), nx_(nx_) {}
    Array1D(const Array1D& other) : detail::ArrayBase<T, 1>(other), nx_(other.nx_) {}
    
    Array1D& operator=(const Array1D& other);
    
    ~Array1D() = default;
    
    // Indexing Function
    [[nodiscard]] static inline size_t get_scid(size_t i) { return i; }

    [[nodiscard]] size_t nx() const { return nx_; }
    static constexpr size_t ny() { return 0u; }
    static constexpr size_t nz() { return 0u; }
    
    // Specialized Accessors
    value_t& operator()(size_t i) { return (*this)[i]; }
    const value_t& operator()(size_t i) const { return (*this)[i]; }

    value_t forward_diff_x(const size_t i) const { return (*this)[get_scid(i + 1)] - (*this)[get_scid(i)]; }
    value_t backward_diff_x(const size_t i) const { return (*this)[get_scid(i)] - (*this)[get_scid(i - 1)]; }

    // Unary Negation
    auto operator-() const {
      Array1D result(nx_);
      for (size_t i = 0; i < this->data.size(); i++) { result[i] = value_t(-1) * (*this)[i]; }
      return result;
    }
    
    // Augmented Assignment
    auto& operator+=(const Array1D& other) {
      assert(nx_ == other.nx_);
      for (size_t i = 0; i < this->data.size(); i++) { (*this)[i] += other.data[i]; }
      return *this;
    }
    
    auto& operator-=(const Array1D& other) {
      assert(nx_ == other.nx_);
      for (size_t i = 0; i < this->data.size(); i++) { (*this)[i] -= other.data[i]; }
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
    const size_t nx_;
    //
  };// end class Array1D
  
  // ----- Array2D -----
  template <typename T>
  struct Array2D : public detail::ArrayBase<T, 2> {
    using value_t = typename detail::ArrayBase<T, 2>::value_t;
    using vector_t = typename detail::ArrayBase<T, 2>::vector_t;
    using dimension_t = typename detail::ArrayBase<T, 2>::dimension_t;
    
    // Constructors & Destructor
    explicit Array2D() = default;
    Array2D(size_t nx_, size_t ny_, value_t fill=0.0) : detail::ArrayBase<T, 2>(nx_ * ny_, fill), nx_(nx_), ny_(ny_) {}
    Array2D(vec2<size_t> dims_, value_t fill=0.0) : Array2D(dims_[0], dims_[1], fill) {}
    Array2D(const Array2D& other) : detail::ArrayBase<T, 2>(other), nx_(other.nx_), ny_(other.ny_) {}
    
    Array2D& operator=(const Array2D& other);
    
    ~Array2D() = default;
    
    // Indexing Function
    [[nodiscard]] inline size_t get_scid(size_t i, size_t k) const { return k + ny_ * i; }

    [[nodiscard]] size_t nx() const { return nx_; }
    [[nodiscard]] size_t ny() const { return ny_; }
    static constexpr size_t nz() { return 0u; }
    
    // Specialized accessors
    value_t& operator()(size_t i, size_t k) { return (*this)[get_scid(i, k)]; }
    const value_t& operator()(size_t i, size_t k) const { return (*this)[get_scid(i, k)]; }

    value_t forward_diff_x(const size_t i, const size_t j) const { return (*this)[get_scid(i + 1, j)] - (*this)[get_scid(i, j)]; }
    value_t backward_diff_x(const size_t i, const size_t j) const { return (*this)[get_scid(i, j)] - (*this)[get_scid(i - 1, j)]; }

    value_t forward_diff_y(const size_t i, const size_t j) const { return (*this)[get_scid(i, j + 1)] - (*this)[get_scid(i, j)]; }
    value_t backward_diff_y(const size_t i, const size_t j) const { return (*this)[get_scid(i, j)] - (*this)[get_scid(i, j - 1)]; }

    // value_t forward_diff_z(const size_t i, const size_t j) const { return (*this)[get_scid(i, j + 1)] - (*this)[get_scid(i, j)]; }
    // value_t backward_diff_z(const size_t i, const size_t j) const { return (*this)[get_scid(i, j)] - (*this)[get_scid(i, j - 1)]; }


    // Dims
    [[nodiscard]] vec2<size_t> dims() const { return {nx_, ny_}; }
    
    // Unary Negation
    auto operator-() const {
      Array2D result(nx_, ny_);
      for (size_t i = 0; i < this->data.size(); i++) { result[i] = value_t(-1) * (*this)[i]; }
      return result;
    }
    
    // Augmented Assignment
    auto& operator+=(const Array2D& other) {
      assert(nx_ == other.nx_ && ny_ == other.ny_);
      for (size_t i = 0; i < this->data.size(); i++) { (*this)[i] += other.data[i]; }
      return *this;
    }
    
    auto& operator-=(const Array2D& other) {
      assert(nx_ == other.nx_ && ny_ == other.ny_);
      for (size_t i = 0; i < this->data.size(); i++) { (*this)[i] -= other.data[i]; }
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
    const size_t nx_, ny_;
    //
  };// end class Array2D
  
  // ----- Array3D -----
  template <typename T>
  struct Array3D : public detail::ArrayBase<T, 3> {
    using value_t = typename detail::ArrayBase<T, 3>::value_t;
    using vector_t = typename detail::ArrayBase<T, 3>::vector_t;
    using dimension_t = typename detail::ArrayBase<T, 3>::dimension_t;
    
    // Constructors & Destructor
    explicit Array3D() = default;
    Array3D(size_t nx_, size_t ny_, size_t nz_, value_t fill=0.0) : detail::ArrayBase<T, 3>(nx_ * ny_ * nz_, fill), nx_(nx_), ny_(ny_), nz_(nz_) {}
    Array3D(vec3<size_t> dims_, value_t fill=0.0) : Array3D(dims_[0], dims_[1], dims_[2], fill) {}
    Array3D(const Array3D& other) : detail::ArrayBase<T, 3>(other), nx_(other.nx_), ny_(other.ny), nz_(other.nz_) {}
    
    Array3D& operator=(const Array3D& other);
    
    ~Array3D() = default;
    
    // Indexing Function
    [[nodiscard]] size_t get_scid(const size_t i, const size_t j, const size_t k) const { return k + nz_ * j + nz_ * ny_ * i; }
    // Dims
    [[nodiscard]] inline tf::types::vec3<size_t> dims() const { return {nx_, ny_, nz_}; }

    [[nodiscard]] size_t nx() const { return nx_; }
    [[nodiscard]] size_t ny() const { return ny_; }
    [[nodiscard]] size_t nz() const { return nz_; }

    // Specialized accessors
    value_t& operator()(size_t i, size_t j, size_t k) { return (*this)[get_scid(i, j, k)]; }
    const value_t& operator()(size_t i, size_t j, size_t k) const { return (*this)[get_scid(i, j, k)]; }
    
    void fill(value_t value) { for (auto& el : this->data) el = value; }

    value_t forward_diff_x(const size_t i, const size_t j, const size_t k) const { return (*this)[get_scid(i + 1, j, k)] - (*this)[get_scid(i, j, k)]; }
    value_t backward_diff_x(const size_t i, const size_t j, const size_t k) const { return (*this)[get_scid(i, j, k)] - (*this)[get_scid(i - 1, j, k)]; }

    value_t forward_diff_y(const size_t i, const size_t j, const size_t k) const { return (*this)[get_scid(i, j + 1, k)] - (*this)[get_scid(i, j, k)]; }
    value_t backward_diff_y(const size_t i, const size_t j, const size_t k) const { return (*this)[get_scid(i, j, k)] - (*this)[get_scid(i, j - 1, k)]; }

    value_t forward_diff_z(const size_t i, const size_t j, const size_t k) const { return (*this)[get_scid(i, j, k + 1)] - (*this)[get_scid(i, j, k)]; }
    value_t backward_diff_z(const size_t i, const size_t j, const size_t k) const { return (*this)[get_scid(i, j, k)] - (*this)[get_scid(i, j, k - 1)]; }

    // Unary Negation
    auto operator-() const {
      Array3D result(nx_, ny_, nz_);
      for (size_t i = 0; i < this->data.size(); i++) { result[i] = value_t(-1) * (*this)[i]; }
      return result;
    }
    
    // Augmented Assignment
    auto& operator+=(const Array3D& other) {
      assert(nx_ == other.nx_ && ny_ == other.ny && nz_ == other.nz_);
      for (size_t i = 0; i < this->data.size(); i++) { (*this)[i] += other.data[i]; }
      return *this;
    }
    
    auto& operator-=(const Array3D& other) {
      assert(nx_ == other.nx_ && ny_ == other.ny && nz_ == other.nz_);
      for (size_t i = 0; i < this->data.size(); i++) { (*this)[i] -= other.data[i]; }
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
    const size_t nx_, ny_, nz_;
    //
  };// end class Array3D
  //
} // end namespace tf::types

// ===== Copy-Assignment Implementations =====
// ===========================================
template <typename T>
tf::types::Array1D<T>& tf::types::Array1D<T>::operator=(const tf::types::Array1D<T> &other) {
  assert(false);
  return *this;
}

template <typename T>
tf::types::Array2D<T>& tf::types::Array2D<T>::operator=(const tf::types::Array2D<T>& other) {
  assert(false);
  return *this;
}

template <typename T>
tf::types::Array3D<T>& tf::types::Array3D<T>::operator=(const tf::types::Array3D<T> &other) {
  assert(false);
  return *this;
}

#endif //TRIFORCE_ARRAY_H
