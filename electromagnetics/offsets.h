//
// Created by cepheid on 10/25/24.
//

#ifndef OFFSETS_H
#define OFFSETS_H

#include <ostream>

struct IntegratorOffsets {
  size_t x0, x1, y0, y1, z0, z1;
};

inline std::ostream& operator<<(std::ostream& os, const IntegratorOffsets& offset) {
  return os << "{" << offset.x0 << ", " << offset.x1 << ", " << offset.y0 << ", " << offset.y1 << ", " << offset.z0 << ", " << offset.z1 << "}";
}

enum class EMFace { X, Y, Z};
enum class EMSide { Lo, Hi};
enum class EMComponent { E, H };

using EMField = EMFace;

#endif //OFFSETS_H
