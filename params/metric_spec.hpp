#ifndef GAME_ENGINE_METRIC_SPEC_HPP
#define GAME_ENGINE_METRIC_SPEC_HPP

#include "mdspan.hpp"

#include <array>
#include <string>

struct FieldSlice {
   std::string component;
   std::array<std::size_t, 3> starts;
   std::array<std::size_t, 3> strides;
   std::extents<std::size_t> extents;
};

#endif //GAME_ENGINE_METRIC_SPEC_HPP