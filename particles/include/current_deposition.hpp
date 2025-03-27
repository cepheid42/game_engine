#ifndef CURRENT_DEPOSITION_HPP
#define CURRENT_DEPOSITION_HPP

#include "math_utils.hpp"

#include "particles.hpp"
#include "em_data.hpp"

#include <array>
#include <cmath>

namespace tf::particles {
  inline std::array<compute_t, 3> quad_shapes(const compute_t v) {
    const auto shape01 = 0.5f * math::SQR(0.5f - v);
    const auto shape2  = 0.75f - math::SQR(v);
    return {shape01, shape01, shape2};
  }


  // inline void split_trajectory(const Particle& p) {
  //   vec3<int> stencil1{static_cast<int>(p.location[0]), static_cast<int>(p.location[1]), static_cast<int>(p.location[2])};
  // }


  struct CurrentDeposition {
    using emdata_t = electromagnetics::EMData;

    static void updateJ(const auto& as0, const auto& as1, const auto& bs0, const auto& bs1, const auto& cs0, const auto& cs1, const auto& qA, const auto& J, const std::size_t i, const std::size_t j, const std::size_t k) {
      static constexpr auto third = 1.0f / 3.0f;
      static constexpr auto sixth = 1.0f / 6.0f;

      auto wT = 0.0f;
      const auto wm = as0[0] - as1[0];
      const auto wp = as1[2] - as0[2];
      for (std::size_t ii = 0; ii < 2; ii++) {
        for (std::size_t jj = 0; jj < 3; jj++) {
          for (std::size_t kk = 0; kk < 3; kk++) {
            wT = qA * (third * (bs0[jj] * cs0[kk] + bs1[jj] * cs1[kk]) + sixth * (bs1[jj] * cs0[kk] + bs0[jj] * cs1[kk]));
            // todo: need to be concerned about edges of grid here, where j - 1 < 0, j + 2 > Ny,... etc
            J(i + ii, j + jj, k + kk) += wm * wT;
          }
        }
      }
    } // end updateJ()

    static void update(const std::vector<Particle>& particles, const compute_t charge, const emdata_t& emdata) {

      for (std::size_t pid = 0; pid < particles.size(); pid++) {
        for (const auto& segment: split_trajectory()) {
          const auto& [i, j, k] = segment.cid;

          const auto xs0 = quad_shapes(segment.p0[0]);
          const auto xs1 = quad_shapes(segment.p1[0]);

          const auto ys0 = quad_shapes(segment.p0[1]);
          const auto ys1 = quad_shapes(segment.p1[1]);

          const auto zs0 = quad_shapes(segment.p0[2]);
          const auto zs1 = quad_shapes(segment.p1[2]);

          const auto cx = segment.weight * charge / (Ayz * dt);
          const auto cy = segment.weight * charge / (Axz * dt);
          const auto cz = segment.weight * charge / (Axy * dt);

          updateJ(xs0, xs1, ys0, ys1, zs0, zs1, cx, emdata.Jx, i, j - 1, k - 1);
          updateJ(ys0, ys1, xs0, xs1, zs0, zs1, cy, emdata.Jy, i - 1, j, k - 1);
          updateJ(zs0, zs1, xs0, xs1, ys0, ys1, cz, emdata.Jz, i - 1, j - 1, k);
        }
      }
    } // end update()
  }; // end struct CurrentDeposition
} // end namepsace tf::particles
#endif //CURRENT_DEPOSITION_HPP
