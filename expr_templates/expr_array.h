//
// Created by cepheid on 8/19/24.
//

#ifndef EXPR_ARRAY_H
#define EXPR_ARRAY_H

#include <cassert>
#include <cstddef>

#include "sarray.h"

template<typename T, typename Rep=SArray<T>>
class Array {
  Rep expr_rep; // (access to) the data of the array

public:
  // create array with initial size
  explicit Array(const std::size_t s) : expr_rep(s) {}
  // explicit Array(const std::size_t sx, const std::size_t sy) : expr_rep(sx, sy) {}

  // create array from possible representation
  explicit Array(const Rep& rb) : expr_rep(rb) {}

  // assignment operator for same type
  Array& operator=(const Array& b) {
    assert(size()==b.size());
    for (std::size_t idx = 0; idx < b.size(); idx++) {
      expr_rep[idx] = b[idx];
    }
    return *this;
  }

  // assignment operator for arrays of different type
  template<typename T2, typename Rep2>
  Array& operator=(const Array<T2, Rep2>& b) {
    assert(size()==b.size());
    for (std::size_t idx = 0; idx < b.size(); idx++) {
      expr_rep[idx] = b[idx];
    }
    return *this;
  }

  // size is size of represented data
  [[nodiscard]] std::size_t size() const { return expr_rep.size(); }

  [[nodiscard]] Shape shape() const { return expr_rep.shape();  }

  // index operator for constants and variables
  decltype(auto) operator[](const std::size_t idx) const {
    assert(idx<size());
    return expr_rep[idx];
  }

  T& operator[](const std::size_t idx) {
    assert(idx<size());
    return expr_rep[idx];
  }

  // // 2D index operator
  // decltype(auto) operator()(const std::size_t i, const std::size_t j) const {
  //   return expr_rep(i, j);
  // }
  //
  // T& operator()(const std::size_t i, const std::size_t j) {
  //   return expr_rep(i, j);
  // }

  template<typename T2, typename R2>
  decltype(auto) operator[](const Array<T2, R2>&) const;

  template<typename T2, typename R2>
  decltype(auto) operator[](const Array<T2, R2>&);

  // template<typename T2, typename R2>
  // decltype(auto) operator()(const Array<T2, R2>&, const Array<T2, R2>&) const;
  //
  // template<typename T2, typename R2>
  // decltype(auto) operator()(const Array<T2, R2>&, const Array<T2, R2>&);

  // return what the array currently represents
  const Rep& rep() const { return expr_rep; }
  Rep& rep() { return expr_rep; }
};

#endif //EXPR_ARRAY_H
