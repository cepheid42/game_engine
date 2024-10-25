//
// Created by cepheid on 10/25/24.
//

#ifndef OFFSETS_H
#define OFFSETS_H

constexpr size_t dPML = 10u;
constexpr size_t nHalo = 2u;

enum class EMFace { X, Y, Z };

struct IntegratorOffsets {
  size_t x0, x1, y0, y1, z0, z1;
};

template<EMFace F, bool=false>
struct PeriodicOffsets {

  explicit PeriodicOffsets(const auto& f) requires (F == EMFace::X)
  : offsets{0, nHalo, 0, f.ny, 0, f.nz}
  {}

  explicit PeriodicOffsets(const auto& f) requires (F == EMFace::Y)
  : offsets{0, f.nx, 0, nHalo, 0, f.nz}
  {}

  explicit PeriodicOffsets(const auto& f) requires (F == EMFace::Z)
  : offsets{0, f.nx, 0, f.ny, 0, nHalo}
  {}

  explicit PeriodicOffsets(auto&) {}

  IntegratorOffsets offsets{};
};

template<EMFace F, bool HI>
struct PMLOffsets {

  explicit PMLOffsets(const auto& f) requires (F == EMFace::X and !HI)
  : offsets{0, dPML, 0, f.ny, 0, f.nz}
  {}

  explicit PMLOffsets(const auto& f) requires (F == EMFace::Y and !HI)
  : offsets{0, f.nx, 0, nHalo, 0, f.nz}
  {}

  explicit PMLOffsets(const auto& f) requires (F == EMFace::Z and !HI)
  : offsets{0, f.nx, 0, f.ny, 0, nHalo}
  {}

  explicit PMLOffsets(const auto& f) requires (F == EMFace::X and HI)
  : offsets{f.nx - dPML, f.nx, 0, f.ny, 0, f.nz}
  {}

  explicit PMLOffsets(const auto& f) requires (F == EMFace::Y and HI)
  : offsets{0, f.nx, 0, nHalo, 0, f.nz}
  {}

  explicit PMLOffsets(const auto& f) requires (F == EMFace::Z and HI)
  : offsets{0, f.nx, 0, f.ny, 0, nHalo}
  {}

  explicit PMLOffsets(auto&) {}

  IntegratorOffsets offsets{};
};

#endif //OFFSETS_H
