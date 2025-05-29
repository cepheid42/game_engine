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
      struct TrajectoryShapeFunc {
         constexpr explicit TrajectoryShapeFunc(StartFunc&& start_func, EndFunc&& end_func)
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

      // void depositJx(auto& J, const auto& p0, const auto& p1, const auto& cids, const auto qA) {
      //    using JitShape = interp::Jit<interp::TSC>;
      //    using CachedShape = interp::Jit<interp::TSC>;
      //    using ReducedShape = interp::Jit<interp::CIC>;
      //    using JitTrajectory = TrajectoryShapeFunc<JitShape, JitShape>;
      //    using ReducedTrajectory = TrajectoryShapeFunc<ReducedShape, ReducedShape>;
      //    using CachedTrajectory = TrajectoryShapeFunc<CachedShape, CachedShape>;
      //
      //    static constexpr auto third = constants::third<compute_t>;
      //
      //    if (p0[0] == p1[0]) { return; } // no deposition in this direction
      //
      //    const auto& [ci, cj, ck] = cids;
      //    const auto& [x0, y0, z0] = p0;
      //    const auto& [x1, y1, z1] = p1;
      //
      //    const JitTrajectory shapeI(JitShape{x0}, JitShape{x1});
      //    const ReducedTrajectory shapeJ(ReducedShape{y0}, ReducedShape{y1});
      //    const CachedTrajectory shapeK(CachedShape{z0}, CachedShape{z1});
      //
      //    for (int j = -1; j <= 0; ++j) {
      //       const auto s0j = shapeJ.S0(j);
      //       const auto dsj = shapeJ.DS(j);
      //
      //       for (int k = -1; k <= 1; ++k) {
      //          const auto s0k = shapeK.S0(k);
      //          const auto dsk = shapeK.DS(k);
      //
      //          const auto tmp = qA * (s0k * s0j + 0.5_fp * (dsk * s0j + s0k * dsj) + third * dsj * dsk);
      //          auto acc = 0.0;
      //          for (int i = -1; i <= 1; ++i) {
      //             acc += shapeI.DS(i) * tmp;
      //             #pragma omp atomic update
      //             J(ci + i, cj + j, ck + k) += acc;
      //          }
      //       }
      //    }
      // }
      //
      // void depositJy(auto& J, const auto& p0, const auto& p1, const auto& cids, const auto qA) {
      //    using JitShape = interp::Jit<interp::TSC>;
      //    using ReducedShape = interp::Jit<interp::CIC>;
      //    using CachedShape = interp::Jit<interp::TSC>;
      //    using JitTrajectory = TrajectoryShapeFunc<JitShape, JitShape>;
      //    using ReducedTrajectory = TrajectoryShapeFunc<ReducedShape, ReducedShape>;
      //    using CachedTrajectory = TrajectoryShapeFunc<CachedShape, CachedShape>;
      //    // (x, y, z) -> (1, 2, 0) -> (y, z, x) | Ex D = 0
      //    // (x, y, z) -> (2, 0, 1) -> (z, x, y) | Ey D = 1
      //    // (x, y, z) ->           -> (x, y, z) | Ez D = 2
      //    static constexpr vec3 b0{-1, 0, -1};
      //    static constexpr vec3 b1{1, 1, 1};
      //    static constexpr auto third = constants::third<compute_t>;
      //
      //    if (p0[1] == p1[1]) {
      //       return;
      //    } // no deposition in this direction
      //
      //    const auto& [ci, cj, ck] = cids;
      //    const auto& [x0, y0, z0] = p0;
      //    const auto& [x1, y1, z1] = p1;
      //
      //    const JitTrajectory shapeI(JitShape{x0}, JitShape{x1});
      //    const ReducedTrajectory shapeJ(ReducedShape{y0}, ReducedShape{y1});
      //    const CachedTrajectory shapeK(CachedShape{z0}, CachedShape{z1});
      //
      //    for (int k = -1; k <= 1; ++k) {
      //       const auto s0k = shapeK.S0(k);
      //       const auto dsk = shapeK.DS(k);
      //
      //       for (int i = -1; i <= 1; ++i) {
      //          const auto s0i = shapeI.S0(i);
      //          const auto dsi = shapeI.DS(i);
      //
      //          const auto tmp = qA * (s0k * s0i + 0.5_fp * (dsk * s0i + s0k * dsi) + third * dsi * dsk);
      //          auto acc = 0.0;
      //          for (int j = -1; j <= -1; ++j) {
      //             const auto W = shapeJ.DS(j) * tmp;
      //             acc += W;
      //             #pragma omp atomic update
      //             J(ci + i, cj + j, ck + k) += acc;
      //          }
      //       }
      //    }
      // }
      //
      // void depositJz(auto& J, const auto& p0, const auto& p1, const auto& cids, const auto qA) {
      //    using JitShape = interp::Jit<interp::TSC>;
      //    using ReducedShape = interp::Jit<interp::CIC>;
      //    using CachedShape = interp::Jit<interp::TSC>;
      //    using JitTrajectory = TrajectoryShapeFunc<JitShape, JitShape>;
      //    using ReducedTrajectory = TrajectoryShapeFunc<ReducedShape, ReducedShape>;
      //    using CachedTrajectory = TrajectoryShapeFunc<CachedShape, CachedShape>;
      //    // (x, y, z) -> (1, 2, 0) -> (y, z, x) | Ex D = 0
      //    // (x, y, z) -> (2, 0, 1) -> (z, x, y) | Ey D = 1
      //    // (x, y, z) ->           -> (x, y, z) | Ez D = 2
      //    static constexpr vec3 b0{-1, 0, -1};
      //    static constexpr vec3 b1{1, 1, 1};
      //    static constexpr auto third = constants::third<compute_t>;
      //
      //    if (p0[2] == p1[2]) {
      //       return;
      //    } // no deposition in this direction
      //
      //    const auto& [ci, cj, ck] = cids;
      //    const auto& [x0, y0, z0] = p0;
      //    const auto& [x1, y1, z1] = p1;
      //
      //    const JitTrajectory shapeI(JitShape{x0}, JitShape{x1});
      //    const ReducedTrajectory shapeJ(ReducedShape{y0}, ReducedShape{y1});
      //    const CachedTrajectory shapeK(CachedShape{z0}, CachedShape{z1});
      //
      //    for (int i = -1; i <= 1; ++i) {
      //       const auto s0i = shapeI.S0(i);
      //       const auto dsi = shapeI.DS(i);
      //
      //       for (int j = -1; j <= 0; ++j) {
      //          const auto s0j = shapeJ.S0(j);
      //          const auto dsj = shapeJ.DS(j);
      //
      //          const auto tmp = qA * (s0j * s0i + 0.5_fp * (dsj * s0i + s0j * dsi) + third * dsi * dsj);
      //          auto acc = 0.0;
      //          for (int k = -1; k <= 1; ++k) {
      //             const auto W = shapeK.DS(k) * tmp;
      //             acc += W;
      //             #pragma omp atomic update
      //             J(ci + i, cj + j, ck + k) += acc;
      //          }
      //       }
      //    }
      // }

      template<int C>
      constexpr auto makeTrajectoryFunction(const auto x0, const auto x1) {
         using CachedTSC = interp::Cached<interp::TSC>;
         using CachedCIC = interp::Cached<interp::CIC>;
         using CachedTrajectory = TrajectoryShapeFunc<CachedTSC, CachedTSC>;
         using ReducedTrajectory = TrajectoryShapeFunc<CachedCIC, CachedCIC>;

         if constexpr (C == 0 or C == 3) {
            return ReducedTrajectory(CachedCIC{x0}, CachedCIC{x1});
         } else {
            return CachedTrajectory(CachedTSC{x0}, CachedTSC{x1});
         }
      }

      template<int D>
      void deposit(auto& J, const auto& p0, const auto& p1, const auto& cids, const auto qA) {
         using CachedTSC = interp::Cached<interp::TSC>;
         using CachedCIC = interp::Cached<interp::CIC>;
         using CachedTrajectory = TrajectoryShapeFunc<CachedTSC, CachedTSC>;
         using ReducedTrajectory = TrajectoryShapeFunc<CachedCIC, CachedCIC>;

         static constexpr auto third = constants::third<compute_t>;
         static constexpr vec3 b0{-1, -1, -1};
         static constexpr vec3 b1 = interp::rotateOrigin<D>(vec3{1, D == 1 ? -1 : 0, 1});

         if (p0[D] == p1[D]) { return; }

         const auto& [ci, cj, ck] = interp::rotateOrigin<D>(cids);
         const auto& [x0, y0, z0] = interp::rotateOrigin<D>(p0);
         const auto& [x1, y1, z1] = interp::rotateOrigin<D>(p1);

         const auto shapeI = makeTrajectoryFunction<D>(x0, x1);
         const auto shapeJ = makeTrajectoryFunction<D + 1>(y0, y1);
         const auto shapeK = makeTrajectoryFunction<D + 2>(z0, z1);

         if constexpr (D == 0) {
            static_assert(std::is_same_v<std::remove_cvref_t<decltype(shapeI)>, ReducedTrajectory>);
            static_assert(std::is_same_v<std::remove_cvref_t<decltype(shapeJ)>, CachedTrajectory>);
            static_assert(std::is_same_v<std::remove_cvref_t<decltype(shapeK)>, CachedTrajectory>);
         }
         else if constexpr (D == 1) {
            static_assert(std::is_same_v<std::remove_cvref_t<decltype(shapeK)>, ReducedTrajectory>);
            static_assert(std::is_same_v<std::remove_cvref_t<decltype(shapeJ)>, CachedTrajectory>);
            static_assert(std::is_same_v<std::remove_cvref_t<decltype(shapeI)>, CachedTrajectory>);
         }
         else {
            static_assert(std::is_same_v<std::remove_cvref_t<decltype(shapeJ)>, ReducedTrajectory>);
            static_assert(std::is_same_v<std::remove_cvref_t<decltype(shapeI)>, CachedTrajectory>);
            static_assert(std::is_same_v<std::remove_cvref_t<decltype(shapeK)>, CachedTrajectory>);
         }

         for (int i = b0[0]; i <= b1[0]; ++i) {
            const auto s0i = shapeI.S0(i);
            const auto dsi = shapeI.S1(i) - s0i;

            for (int j = b0[1]; j <= b1[1]; ++j) {
               const auto s0j = shapeJ.S0(j);
               const auto dsj = shapeJ.S1(j) - s0j;

               const auto tmp = qA * (s0i * s0j + 0.5_fp * (dsi * s0j + s0i * dsj) + third * dsj * dsi);
               auto acc = 0.0;
               for (int k = b0[2]; k <= b1[2]; ++k) {
                  acc += shapeK.DS(k) * tmp;

                  const auto [x, y, z] = interp::rotateOrigin<D == 2 ? D : !D>(vec3{ci + i, cj + j, ck + k});
                  #pragma omp atomic update
                  J(x, y, z) += acc;
               }
            }
         }
      }

      auto findRelayPoint(auto& i0, auto& i1, const auto& x0, const auto& x1, const bool is_even) {
         if (is_even) {
            i0 = std::floor(x0);
            i1 = std::floor(x1);
            return i0 == i1 ? x1 : std::max(i0, i1);
         }

         i0 = std::floor(x0 + 0.5_fp);
         i1 = std::floor(x1 + 0.5_fp);
         return i0 == i1 ? x1 : (i0 + i1) / 2.0_fp;
      }

      void updateJ(const auto& p, auto& emdata, const auto charge) {
         if (p.disabled) { return; }

         const auto x_coeff = p.weight * charge * dtAyz;
         const auto z_coeff = p.weight * charge * dtAxy;
         const auto y_coeff = p.weight * charge * dtAxz;

         auto cids = getCIDs(p.old_location + 0.5_fp);

         // const auto pold = p.old_location;
         // const auto pnew = p.location;

         vec3<compute_t> i0{}, i1{}, relay{};
         for (int d = 0; d < 3; ++d) {
            // p0[d] = pold[d];
            relay[d] = findRelayPoint(i0[d], i1[d], p.old_location[d], p.location[d], d == 1);
         }

         auto p0 = p.old_location - i0;
         auto p1 = relay - i0;

         deposit<0>(emdata.Jx, p0, p1, cids, x_coeff);
         deposit<1>(emdata.Jy, p0, p1, cids, y_coeff);
         deposit<2>(emdata.Jz, p0, p1, cids, z_coeff);

         if (i0 != i1) {
            for (int d = 0; d < 3; ++d) {
               p0[d] = relay[d] - i1[d];
               p1[d] = p.location[d] - i1[d];
               cids[d] = cids[d] + (i1[d] - i0[d]);
            }

            deposit<0>(emdata.Jx, p0, p1, cids, x_coeff);
            // deposit<1>(emdata.Jy, p0, p1, cids, y_coeff);
            deposit<2>(emdata.Jz, p0, p1, cids, z_coeff);
         }
      }
   } // namespace detail

   void deposit_current(auto& emdata, const auto& g) {
      #pragma omp parallel for num_threads(nThreads)
      for (std::size_t pid = 0; pid < g.num_particles(); pid++) {
         detail::updateJ(g.particles[pid], emdata, g.charge);
      }
   }
} // namespace tf::particles

#endif  // CURRENT_DEPOSITION_HPP
