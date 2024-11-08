//
// Created by cepheid on 10/21/24.
//

#ifndef CURL_OPERATORS_H
#define CURL_OPERATORS_H

#include "electromagnetics.param"


//====== Curl Operators =======
//=============================
template<Derivative D, bool Forward, typename... IDXS>
struct curl {
  static constexpr auto apply(const auto&, IDXS...) {
    return 0.0;
  }
};

template<bool Forward, typename... IDXS>
struct curl<Derivative::DX, Forward, IDXS...> {
  static auto apply(const auto& f, IDXS... idxs) {
    if constexpr (Forward) {
      return f.forward_diff_x(idxs...);
    } else {
      return f.backward_diff_x(idxs...);
    }
  }
};

template<bool Forward, typename... IDXS>
struct curl<Derivative::DY, Forward, IDXS...> {
  static auto apply(const auto& f, IDXS... idxs) {
    if constexpr (Forward) {
      return f.forward_diff_y(idxs...);
    } else {
      return f.backward_diff_y(idxs...);
    }
  }
};

template<bool Forward, typename... IDXS>
struct curl<Derivative::DZ, Forward, IDXS...> {
  static auto apply(const auto& f, IDXS... idxs) {
    if constexpr (Forward) {
      return f.forward_diff_z(idxs...);
    } else {
      return f.backward_diff_z(idxs...);
    }
  }
};

#endif //CURL_OPERATORS_H
