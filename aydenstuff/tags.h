//
// Created by akis on 8/27/24.
//


#ifndef TRIFORCE_TAGS_H
#define TRIFORCE_TAGS_H

#include <concepts>

using std::size_t;

// ===== Tags =====
// ================
namespace tf::tags
{
  template <size_t N>
  struct Dimension {
    static constexpr size_t value = N;
  };
  //
} // end namespace tf::tags

namespace tf::concepts
{
  template <typename T>
  concept is_dimensional = requires { typename T::dimension_t; };
  
  template <typename T, size_t N>
  concept is_n_dimensional = std::same_as<typename T::dimension_t, tags::Dimension<N>>;

  template<typename T1, typename T2>
  concept is_same_dimensional = std::same_as<typename T1::dimension_t, typename T2::dimension_t>;
  //
} // end namespace tf::concepts

#endif //TRIFORCE_TAGS_H
