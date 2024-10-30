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
  // ----- Dimension -----
  template <size_t N>
  struct Dimension { static constexpr size_t value = N; };
  
  // ----- Order -----
  template <size_t N>
  struct Order { static constexpr size_t value = N; };
  //
} // end namespace tf::tags

namespace tf::concepts
{
  // ----- Dimension -----
  template <typename T>
  concept has_dimension = requires { typename T::dimension_t; };
  
  template <typename T, size_t N>
  concept is_dimension_n = std::same_as<typename T::dimension_t, tags::Dimension<N>>;
  
  // ----- Order -----
  template <typename T>
  concept has_order = requires { typename T::order_t; };
  
  template <typename T, size_t N>
  concept is_order_n = std::same_as<typename T::dimension_t, tags::Order<N>>;
  //
} // end namespace tf::concepts

#endif //TRIFORCE_TAGS_H
