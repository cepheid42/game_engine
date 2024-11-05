//
// Created by cepheid on 10/23/24.
//

#ifndef BC_DATA_H
#define BC_DATA_H

// #include "em_traits.h"
// #include "aydenstuff/array.h"
#include "boundaries.h"

using tf::types::Array1D;
using tf::types::Array2D;
using tf::types::Array3D;

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

template<typename EzX0, typename EzX1>
struct BCData {

  explicit BCData(const auto& Ez)
  : ez_x0{get_offsets<EMFace::X, EMSide::Lo, EzX0::bc_depth>(Ez)},
    ez_x1{get_offsets<EMFace::X, EMSide::Hi, EzX1::bc_depth>(Ez)}
  {}

  EzX0 ez_x0;
  EzX1 ez_x1;
};


#endif //BC_DATA_H
