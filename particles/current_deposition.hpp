#ifndef CURRENT_DEPOSITION_HPP
#define CURRENT_DEPOSITION_HPP

#include <cmath>

#include "constants.hpp"
#include "em_data.hpp"
#include "interpolation.hpp"
#include "particles.hpp"

namespace tf::particles {
namespace detail {
template <typename StartFunc, typename EndFunc>
struct TrajectoryAssignment {
  constexpr explicit TrajectoryAssignment(StartFunc&& start_func,
                                          EndFunc&& end_func)
      : startShape(std::move(start_func)), endShape(std::move(end_func)) {}

  [[nodiscard]] constexpr auto S0(const int grid_point) const {
    return startShape(grid_point);
  }

  [[nodiscard]] constexpr auto S1(const int grid_point) const {
    return endShape(grid_point);
  }

  [[nodiscard]] constexpr auto DS(const int grid_point) const {
    return S1(grid_point) - S0(grid_point);
  }

  const StartFunc startShape;
  const EndFunc endShape;
};

void depositJx(auto& J, const auto& p0, const auto& p1, const auto& cids, const auto qA) {
  using JitShape = interp::Jit<interp::TSC>;
  using CachedShape = interp::Jit<interp::TSC>;
  using ReducedShape = interp::Jit<interp::CIC>;
  using JitTrajectory = TrajectoryAssignment<JitShape, JitShape>;
  using ReducedTrajectory = TrajectoryAssignment<ReducedShape, ReducedShape>;
  using CachedTrajectory = TrajectoryAssignment<CachedShape, CachedShape>;
  // (x, y, z) -> (1, 2, 0) -> (y, z, x) | Ex D = 0
  // (x, y, z) -> (2, 0, 1) -> (z, x, y) | Ey D = 1
  // (x, y, z) ->           -> (x, y, z) | Ez D = 2
  static constexpr vec3 b0{-1, 0, -1};
  static constexpr vec3 b1{1, 1, 1};
  static constexpr auto third = constants::third<compute_t>;

  if (p0[0] == p1[0]) {
    return;
  }  // no deposition in this direction

  const auto& [ci, cj, ck] = cids;
  const auto& [x0, y0, z0] = p0;
  const auto& [x1, y1, z1] = p1;

  const JitTrajectory shapeI(JitShape{x0}, JitShape{x1});
  const ReducedTrajectory shapeJ(ReducedShape{y0}, ReducedShape{y1});
  const CachedTrajectory shapeK(CachedShape{z0}, CachedShape{z1});

  for (int j = -1; j <= 0; ++j) {
    const auto s0j = shapeJ.S0(j);
    const auto dsj = shapeJ.DS(j);

    for (int k = -1; k <= 1; ++k) {
      const auto s0k = shapeK.S0(k);
      const auto dsk = shapeK.DS(k);

      const auto tmp = qA * (s0k * s0j + 0.5_fp * (dsk * s0j + s0k * dsj) + third * dsj * dsk);
      auto acc = 0.0;
      for (int i = -1; i <= 1; ++i) {
        const auto W = shapeI.DS(i) * tmp;
        acc += W;
#pragma omp atomic update
        J(ci + i, cj + j, ck + k) += acc;
      }
    }
  }
}

void depositJy(auto& J, const auto& p0, const auto& p1, const auto& cids, const auto qA) {
  using JitShape = interp::Jit<interp::TSC>;
  using ReducedShape = interp::Jit<interp::CIC>;
  using CachedShape = interp::Jit<interp::TSC>;
  using JitTrajectory = TrajectoryAssignment<JitShape, JitShape>;
  using ReducedTrajectory = TrajectoryAssignment<ReducedShape, ReducedShape>;
  using CachedTrajectory = TrajectoryAssignment<CachedShape, CachedShape>;
  // (x, y, z) -> (1, 2, 0) -> (y, z, x) | Ex D = 0
  // (x, y, z) -> (2, 0, 1) -> (z, x, y) | Ey D = 1
  // (x, y, z) ->           -> (x, y, z) | Ez D = 2
  static constexpr vec3 b0{-1, 0, -1};
  static constexpr vec3 b1{1, 1, 1};
  static constexpr auto third = constants::third<compute_t>;

  if (p0[1] == p1[1]) {
    return;
  }  // no deposition in this direction

  const auto& [ci, cj, ck] = cids;
  const auto& [x0, y0, z0] = p0;
  const auto& [x1, y1, z1] = p1;

  const JitTrajectory shapeI(JitShape{x0}, JitShape{x1});
  const ReducedTrajectory shapeJ(ReducedShape{y0}, ReducedShape{y1});
  const CachedTrajectory shapeK(CachedShape{z0}, CachedShape{z1});

  for (int k = -1; k <= 1; ++k) {
    const auto s0k = shapeK.S0(k);
    const auto dsk = shapeK.DS(k);

    for (int i = -1; i <= 1; ++i) {
      const auto s0i = shapeI.S0(i);
      const auto dsi = shapeI.DS(i);

      const auto tmp = qA * (s0k * s0i + 0.5_fp * (dsk * s0i + s0k * dsi) + third * dsi * dsk);
      auto acc = 0.0;
      for (int j = -1; j <= -1; ++j) {
        const auto W = shapeJ.DS(j) * tmp;
        acc += W;
#pragma omp atomic update
        J(ci + i, cj + j, ck + k) += acc;
      }
    }
  }
}

  void depositJz(auto& J, const auto& p0, const auto& p1, const auto& cids, const auto qA) {
  using JitShape = interp::Jit<interp::TSC>;
  using ReducedShape = interp::Jit<interp::CIC>;
  using CachedShape = interp::Jit<interp::TSC>;
  using JitTrajectory = TrajectoryAssignment<JitShape, JitShape>;
  using ReducedTrajectory = TrajectoryAssignment<ReducedShape, ReducedShape>;
  using CachedTrajectory = TrajectoryAssignment<CachedShape, CachedShape>;
  // (x, y, z) -> (1, 2, 0) -> (y, z, x) | Ex D = 0
  // (x, y, z) -> (2, 0, 1) -> (z, x, y) | Ey D = 1
  // (x, y, z) ->           -> (x, y, z) | Ez D = 2
  static constexpr vec3 b0{-1, 0, -1};
  static constexpr vec3 b1{1, 1, 1};
  static constexpr auto third = constants::third<compute_t>;

  if (p0[2] == p1[2]) {
    return;
  }  // no deposition in this direction

  const auto& [ci, cj, ck] = cids;
  const auto& [x0, y0, z0] = p0;
  const auto& [x1, y1, z1] = p1;

  const JitTrajectory shapeI(JitShape{x0}, JitShape{x1});
  const ReducedTrajectory shapeJ(ReducedShape{y0}, ReducedShape{y1});
  const CachedTrajectory shapeK(CachedShape{z0}, CachedShape{z1});

  for (int i = -1; i <= 1; ++i) {
    const auto s0i = shapeI.S0(i);
    const auto dsi = shapeI.DS(i);

    for (int j = -1; j <= 0; ++j) {
      const auto s0j = shapeJ.S0(j);
      const auto dsj = shapeJ.DS(j);

      const auto tmp = qA * (s0j * s0i + 0.5_fp * (dsj * s0i + s0j * dsi) + third * dsi * dsj);
      auto acc = 0.0;
      for (int k = -1; k <= 1; ++k) {
        const auto W = shapeK.DS(k) * tmp;
        acc += W;
#pragma omp atomic update
        J(ci + i, cj + j, ck + k) += acc;
      }
    }
  }
}

// void depositJy(auto& Jy, const auto& p0, const auto& p1, const auto& cids,
// const auto qA)
// {
//    using CachedShape           = interp::Jit<interp::TSC>;
//    using CachedTrajectory      = TrajectoryAssignment<CachedShape,
//    CachedShape>; static constexpr auto third = constants::third<compute_t>;
//
//    if (p0[1] == p1[1]) { return; } // No deposition in this direction
//
//    const CachedTrajectory shapeI(CachedShape{p0[0]}, CachedShape{p1[0]});
//    const CachedTrajectory shapeK(CachedShape{p0[2]}, CachedShape{p1[2]});
//
//    const auto& [ci, cj, ck] = cids;
//    for (int i = -1; i <= 1; i++)
//    {
//       const auto s0i = shapeI.S0(ci + i);
//       const auto dsi = shapeI.DS(ci + i);
//       for (int k = -1; k <= 1; k++)
//       {
//          const auto s0k = shapeK.S0(ck + k);
//          const auto dsk = shapeK.DS(ck + k);
//
//          // Does this need a minus sign like the normal version?
//          #pragma omp atomic update
//          Jy(ci + i, cj, ck + k) += qA * third * (s0i * s0k + 0.5_fp * (dsi *
//          s0k + s0i * dsk) + third * dsk * dsi);
//       }
//    }
// }

auto findRelayPoint(auto& i0, auto& i1, const auto& x0, const auto& x1,
                    const bool is_even) {
  if (is_even) {
    i0 = std::floor(x0);
    i1 = std::floor(x1);
    return i0 == i1 ? x1 : std::max(i0, i1);
  }

  i0 = static_cast<int>(std::floor(x0 + 0.5_fp));
  i1 = static_cast<int>(std::floor(x1 + 0.5_fp));
  return i0 == i1 ? x1 : static_cast<compute_t>(i0 + i1) / 2.0_fp;
}

vec3<compute_t> getCellNormLoc(const auto& loc) {
  return {loc[0] - std::floor(loc[0]), loc[1] - std::floor(loc[1]), loc[2] - std::floor(loc[2])};
}

void updateJ(const auto& p, auto& emdata, const auto charge) {
  if (p.disabled) {
    return;
  }

  const auto x_coeff = p.weight * charge / (Ayz * dt);
  const auto z_coeff = p.weight * charge / (Axy * dt);
  const auto y_coeff = p.weight * charge / (Axz * dt);

  auto cids = getCIDs(p.old_location + 0.5);

  const auto pold = p.old_location;
  const auto pnew = p.location;

  vec3<compute_t> i0{}, i1{}, p0{}, p1{}, relay{};
  for (int d = 0; d < 3; ++d) {
    p0[d] = pold[d];
    relay[d] = findRelayPoint(i0[d], i1[d], pold[d], pnew[d], d == 1);
  }

  p0 = p0 - i0;
  p1 = relay - i0;

  depositJx(emdata.Jx, p0, p1, cids, x_coeff);
  depositJy(emdata.Jy, p0, p1, cids, y_coeff);
  depositJz(emdata.Jz, p0, p1, cids, z_coeff);

  if (i0 != i1) {
    for (int d = 0; d < 3; ++d) {
      p0[d] = relay[d] - i1[d];
      p1[d] = pnew[d] - i1[d];
      cids[d] = cids[d] + (i1[d] - i0[d]);
    }

    depositJx(emdata.Jx, p0, p1, cids, x_coeff);
    // depositJy(emdata.Jy, p0, p1, cids, y_coeff);
    depositJz(emdata.Jz, p0, p1, cids, z_coeff);
  }
}
}  // namespace detail

void deposit_current(auto& emdata, const auto& g) {
  // #pragma omp parallel for num_threads(nThreads)
  for (std::size_t pid = 0; pid < g.num_particles(); pid++) {
    detail::updateJ(g.particles[pid], emdata, g.charge);
  }
}
}  // namespace tf::particles

#endif  // CURRENT_DEPOSITION_HPP
