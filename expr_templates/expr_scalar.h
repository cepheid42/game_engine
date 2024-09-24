//
// Created by cepheid on 8/19/24.
//

#ifndef EXPR_SCALAR_H
#define EXPR_SCALAR_H

#include <cassert>
#include <cstddef>

// class for objects that represent scalars:
template<typename T>
class A_Scalar {
  const T& s; // value of the scalar

public:
  // constructor initializes value
  constexpr explicit A_Scalar(const T& v) : s(v) {}

  // for index operations, the scalar is the value of each element
  constexpr const T& operator[](std::size_t) const { return s; }

  // scalars have zero as size
  [[nodiscard]] static constexpr std::size_t size() { return 0; }
};


#endif //EXPR_SCALAR_H
