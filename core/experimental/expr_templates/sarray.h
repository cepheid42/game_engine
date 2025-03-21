//
// Created by cepheid on 8/19/24.
//

#ifndef SARRAY_H
#define SARRAY_H

#include <cassert>
#include <cstddef>
#include <memory>

struct Shape {
  std::size_t x, y, z;
};

template<typename T>
class SArray {
  std::unique_ptr<T[]> m_data;
  std::size_t m_size;
  Shape m_shape;

public:
  // create array with initial size
  explicit SArray(const std::size_t s)
  : m_data(new T[s]{}),
    m_size(s),
    m_shape{s, 0, 0}
  {}

  // explicit SArray(const std::size_t sx, const std::size_t sy)
  // : m_data(new T[sx * sy]{}),
  //   m_size(sx * sy),
  //   m_shape{sx, sy, 0}
  // {}

  // explicit SArray(const std::size_t sx, const std::size_t sy, const std::size_t sz)
  // : m_data(new T[sx * sy * sz]{}),
  //   m_size(sx * sy * sz),
  //   m_shape{sx, sy, sz}
  // {}

  // copy constructor
  SArray(const SArray& orig)
  : m_data(new T[orig.size()]),
    m_size(orig.size()),
    m_shape(orig.shape())
  { copy(orig); }

  // destructor: free memory
  ~SArray() = default;

  // assignment operator
  SArray& operator=(const SArray& orig) {
    if (&orig != this) {
      copy(orig);
    }
    return *this;
  }

  // return size
  [[nodiscard]] std::size_t size() const { return m_size; }
  // return shape
  [[nodiscard]] Shape shape() const { return m_shape; }

  // [[nodiscard]] std::size_t index(const std::size_t i, const std::size_t j) const { return j + (i * m_shape.y); }

  // index operator for constants and variables
  T& operator[](const std::size_t idx) { return m_data[idx]; }
  const T& operator[](const std::size_t idx) const { return m_data[idx]; }

  // T& operator()(const std::size_t i, const std::size_t j) { return m_data[index(i, j)]; }
  // const T& operator()(const std::size_t i, const std::size_t j) const { return m_data[index(i, j)]; }

protected:
  // copy values of another array
  void copy(const SArray& orig) {
    assert(size() == orig.size());
#pragma omp parallel for
    for (std::size_t idx = 0; idx < size(); ++idx) {
      m_data[idx] = orig.m_data[idx];
    }
  }
};


#endif //SARRAY_H
