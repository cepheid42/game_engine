//
// Created by cepheid on 10/17/24.
//

#ifndef EM_EMTPYARRAY_H
#define EM_EMTPYARRAY_H

template<typename T, std::size_t N>
struct EmptyArray {
  using value_t = T;
  // using vector_t = std::vector<value_t>;
  using dimension_t = tf::tags::Dimension<N>;

  EmptyArray() = default;
  explicit EmptyArray(std::size_t...) {}

  constexpr value_t operator[](std::size_t) const { return static_cast<value_t>(0.0); }
  constexpr value_t operator()(std::size_t...) const { return static_cast<value_t>(0.0); }
};

template<typename T>
using EmptyArray1D = EmptyArray<T, 1>;

template<typename T>
using EmptyArray2D = EmptyArray<T, 2>;

template<typename T>
using EmptyArray3D = EmptyArray<T, 3>;


#endif //EM_EMTPYARRAY_H
