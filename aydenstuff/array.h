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
    // =================== Empty Array Class for Electromagnetics =======================
    // ==================================================================================
    template<typename Array>
    struct EmptyArray {
      using value_t = typename Array::value_t;
      using dimension_t = typename Array::dimension_t;
      using array_t = Array;

      EmptyArray() = default;
      explicit EmptyArray(std::size_t...) {}

      static constexpr size_t nx() { return 0ul; }
      static constexpr size_t ny() { return 0ul; }
      static constexpr size_t nz() { return 0ul; }

      constexpr value_t operator[](std::size_t) const { return static_cast<value_t>(0.0); }
      constexpr value_t operator()(std::size_t...) const { return static_cast<value_t>(0.0); }
    };

    template <typename T, size_t N>
    class ArrayBase {
    public:
      using value_t = T;
      using vector_t = std::vector<value_t>;
      using dimension_t = tags::Dimension<N>;
      
      explicit ArrayBase() = default;
      explicit ArrayBase(const size_t size) : data(size) {}
      ArrayBase(const size_t size, const value_t fill) : data(size, fill) {}
      ArrayBase(const ArrayBase& other) : data(other.data) {}
      
      ~ArrayBase() = default;
      
      [[nodiscard]] size_t size() const { return data.size(); }
      [[nodiscard]] size_t num_bytes() const { return size() * sizeof(value_t); }
      
      explicit operator bool() const { return data != nullptr; }

      auto begin() { return data.begin(); }
      auto end() { return data.end(); }
      
      void fill(value_t value) { for (auto& el : this->data) el = value; }

      // Basic accessors
      value_t& operator[](const size_t i) { return data[i]; }
      const value_t& operator[](const size_t i) const { return data[i]; }
    
    protected:
      vector_t data;
    };
  }
  
  // ----- Array1D -----
  template <typename T>
  struct Array1D : detail::ArrayBase<T, 1> {
    using value_t = typename detail::ArrayBase<T, 1>::value_t;
    using vector_t = typename detail::ArrayBase<T, 1>::vector_t;
    using dimension_t = typename detail::ArrayBase<T, 1>::dimension_t;
    using array_t = Array1D;

    // Constructors & Destructor
    explicit Array1D() = default;
    explicit Array1D(const size_t nx_) : detail::ArrayBase<T, 1>(nx_), nx_(nx_) {}
    Array1D(const size_t nx_, const value_t fill) : detail::ArrayBase<T, 1>(nx_, fill), nx_(nx_) {}
    Array1D(const Array1D& other) : detail::ArrayBase<T, 1>(other), nx_(other.nx()) {}

    Array1D& operator=(const Array1D& other);

    ~Array1D() = default;

    // Indexing Function
    [[nodiscard]] static size_t get_scid(const size_t i) { return i; }

    // Specialized Accessors
    value_t& operator()(size_t i) { return (*this)[i]; }
    const value_t& operator()(size_t i) const { return (*this)[i]; }

    [[nodiscard]] size_t nx() const { return nx_; }

    // Unary Negation
    auto operator-() const {
      Array1D result(nx_);
      for (size_t i = 0; i < this->data.size(); i++) { result[i] = value_t(-1) * (*this)[i]; }
      return result;
    }

    // Augmented Assignment
    auto& operator+=(const Array1D& other) {
      assert(nx_ == other.nx);
      for (size_t i = 0; i < this->data.size(); i++) { (*this)[i] += other.data[i]; }
      return *this;
    }

    auto& operator-=(const Array1D& other) {
      assert(nx_ == other.nx);
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
    size_t nx_;
    //
  };// end class Array1D

  // ----- Array2D -----
  template <typename T>
  struct Array2D : public detail::ArrayBase<T, 2> {
    using value_t = typename detail::ArrayBase<T, 2>::value_t;
    using vector_t = typename detail::ArrayBase<T, 2>::vector_t;
    using dimension_t = typename detail::ArrayBase<T, 2>::dimension_t;
    using array_t = Array2D;

    // Constructors & Destructor
    explicit Array2D() : detail::ArrayBase<T, 2>(), nx_{0ul}, nz_{0ul} {}
    Array2D(const size_t nx_, const size_t nz_, const value_t fill=0.0) : detail::ArrayBase<T, 2>(nx_ * nz_, fill), nx_(nx_), nz_(nz_) {}
    explicit Array2D(const vec2<size_t> dims_, const value_t fill=0.0) : Array2D(dims_[0], dims_[1], fill) {}
    Array2D(const Array2D& other) : detail::ArrayBase<T, 2>(other), nx_(other.nx()), nz_(other.nz()) {}

    Array2D& operator=(const Array2D& other);

    ~Array2D() = default;

    // Indexing Function
    [[nodiscard]] size_t get_scid(const size_t i, const size_t k) const { return k + nz_ * i; }

    // Specialized accessors
    value_t& operator()(const size_t i, const size_t k) { return (*this)[get_scid(i, k)]; }
    const value_t& operator()(const size_t i, const size_t k) const { return (*this)[get_scid(i, k)]; }

    // Dims
    [[nodiscard]] vec2<size_t> dims() const { return {nx_, nz_}; }

    [[nodiscard]] size_t nx() const { return nx_; }
    [[nodiscard]] size_t nz() const { return nz_; }

    // Unary Negation
    auto operator-() const {
      Array2D result(nx_, nz_);
      for (size_t i = 0; i < this->data.size(); i++) { result[i] = value_t(-1) * (*this)[i]; }
      return result;
    }

    // Augmented Assignment
    auto& operator+=(const Array2D& other) {
      assert(nx_ == other.nx && nz_ == other.nz);
      for (size_t i = 0; i < this->data.size(); i++) { (*this)[i] += other.data[i]; }
      return *this;
    }

    auto& operator-=(const Array2D& other) {
      assert(nx_ == other.nx && nz_ == other.nz);
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
    size_t nx_, nz_;
    //
  };// end class Array2D
  
  // ----- Array3D -----
  template <typename T>
  struct Array3D : public detail::ArrayBase<T, 3> {
    using value_t = typename detail::ArrayBase<T, 3>::value_t;
    using vector_t = typename detail::ArrayBase<T, 3>::vector_t;
    using dimension_t = typename detail::ArrayBase<T, 3>::dimension_t;
    using array_t = Array3D;

    // Constructors & Destructor
    explicit Array3D() : detail::ArrayBase<T, 3>(), nx_(0ul), ny_(0ul), nz_(0ul) {}
    explicit Array3D(const size_t nx_, const size_t ny_, const size_t nz_, const value_t fill=0.0) : detail::ArrayBase<T, 3>(nx_ * ny_ * nz_, fill), nx_(nx_), ny_(ny_), nz_(nz_) {}
    // explicit Array3D(vec3<size_t> dims_, value_t fill=0.0) : Array3D(dims_[0], dims_[1], dims_[2], fill) {}
    Array3D(const Array3D& other) : detail::ArrayBase<T, 3>(other), nx_(other.nx()), ny_(other.ny()), nz_(other.nz()) {}
    
    Array3D& operator=(const Array3D& other);
    
    ~Array3D() = default;
    
    // Indexing Function
    [[nodiscard]] size_t get_scid(const size_t i, const size_t j, const size_t k) const { return k + nz_ * j + nz_ * ny_ * i; }

    // Specialized accessors
    value_t& operator()(const size_t i, const size_t j, const size_t k) { return (*this)[get_scid(i, j, k)]; }
    const value_t& operator()(const size_t i, const size_t j, const size_t k) const { return (*this)[get_scid(i, j, k)]; }
    
    void fill(value_t value) { for (auto& el : this->data) el = value; }

    // Dims
    [[nodiscard]] tf::types::vec3<size_t> dims() const { return {nx_, ny_, nz_}; }

    [[nodiscard]] size_t nx() const { return nx_; }
    [[nodiscard]] size_t ny() const { return ny_; }
    [[nodiscard]] size_t nz() const { return nz_; }

    // Unary Negation
    auto operator-() const {
      Array3D result(nx_, ny_, nz_);
      for (size_t i = 0; i < this->data.size(); i++) { result[i] = value_t(-1) * (*this)[i]; }
      return result;
    }
    
    // Augmented Assignment
    auto& operator+=(const Array3D& other) {
      assert(nx_ == other.nx && ny_ == other.ny && nz_ == other.nz);
      for (size_t i = 0; i < this->data.size(); i++) { (*this)[i] += other.data[i]; }
      return *this;
    }
    
    auto& operator-=(const Array3D& other) {
      assert(nx_ == other.nx && ny_ == other.ny && nz_ == other.nz);
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
    size_t nx_, ny_, nz_;
    //
  };// end class Array3D

  template<typename T>
  using EmptyArray1D = detail::EmptyArray<Array1D<T>>;

  template<typename T>
  using EmptyArray2D = detail::EmptyArray<Array2D<T>>;

  template<typename T>
  using EmptyArray3D = detail::EmptyArray<Array3D<T>>;

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
