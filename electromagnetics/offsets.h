//
// Created by cepheid on 10/25/24.
//

#ifndef OFFSETS_H
#define OFFSETS_H

#include "em_traits.h"

namespace tf::electromagnetics::types
{
  using namespace tf::electromagnetics::traits;

  struct IntegratorOffsets {
    size_t x0, x1, y0, y1, z0, z1;
  };

  template<EMSide Side, size_t DEPTH>
  IntegratorOffsets get_x_offsets(const auto& f) {
    if constexpr (Side == EMSide::Lo) {
      return {0, DEPTH, 0, f.ny(), 0, f.nz()};
    } else {
      return {f.nx() - DEPTH, f.nx(), 0, f.ny(), 0, f.nz()};
    }
  }

  template<EMSide Side, size_t DEPTH>
  IntegratorOffsets get_y_offsets(const auto& f) {
    if constexpr (Side == EMSide::Lo) {
      return {0, f.nx(), 0, DEPTH, 0, f.nz()};
    } else {
      return {0, f.nx(), f.ny() - DEPTH, f.ny(), 0, f.nz()};
    }
  }

  template<EMSide Side, size_t DEPTH>
  IntegratorOffsets get_z_offsets(const auto& f) {
    if constexpr (Side == EMSide::Lo) {
      return {0, f.nx(), 0, f.ny(), 0, DEPTH};
    } else {
      return {0, f.nx(), 0, f.ny(), f.nz() - DEPTH, f.nz()};
    }
  }

  template<EMFace Face, EMSide Side, size_t DEPTH>
  IntegratorOffsets get_offsets(const auto& f) {
    if constexpr (Face == EMFace::X) {
      return get_x_offsets<Side, DEPTH>(f);
    } else if constexpr (Face == EMFace::Y) {
      return get_y_offsets<Side, DEPTH>(f);
    } else {
      return get_z_offsets<Side, DEPTH>(f);
    }
  }
}
#endif //OFFSETS_H
