//
// Created by cepheid on 10/25/24.
//

#ifndef OFFSETS_H
#define OFFSETS_H

constexpr size_t dPML = 0u;
constexpr size_t nHalo = 2u;

struct IntegratorOffsets {
  size_t x0, x1, y0, y1, z0, z1;
};

#endif //OFFSETS_H
