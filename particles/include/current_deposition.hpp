#ifndef CURRENT_DEPOSITION_HPP
#define CURRENT_DEPOSITION_HPP

#include "math_utils.hpp"
#include "particles.hpp"
#include "em_data.hpp"
#include "bc_data.hpp"

#include <array>
#include <cmath>

namespace tf::particles {
  inline std::array<compute_t, 3> quad_shapes(const compute_t v) {
    return {
      0.5_fp * math::SQR(0.5_fp - v),
      0.75_fp - math::SQR(v),
      0.5_fp * math::SQR(0.5_fp + v)
    };
  }

  struct Segment {
    std::array<std::size_t, 3> cids = {0, 0, 0};
    vec3<compute_t> p0{};
    vec3<compute_t> p1{};
    bool active{false};
  };

  inline std::array<Segment, 2> split_trajectory(const Particle& p, const std::array<std::size_t, 3>& cidx1) {
    const auto xoff = std::floor(p.old_location[0]);
    const auto yoff = std::floor(p.old_location[1]);
    const auto zoff = std::floor(p.old_location[2]);

    if (xoff == 0.0 and yoff == 0.0 and zoff == 0.0) {
      // Particle did not leave current cell
      return {Segment{cidx1, p.old_location, p.location, true},
              Segment{cidx1, p.old_location, p.location, false}};
    }

    const auto xr0 = xoff == 0.0_fp ? p.location[0] : std::max(0.0_fp, -xoff);
    const auto yr0 = yoff == 0.0_fp ? p.location[1] : std::max(0.0_fp, -yoff);
    const auto zr0 = zoff == 0.0_fp ? p.location[2] : std::max(0.0_fp, -zoff);

    const vec3 pr0{xr0, yr0, zr0};
    const vec3 pr1{xr0 + xoff, yr0 + yoff, zr0 + zoff};

    // todo: pray to the gods that none of the cell id's are zero...
    const auto xo = static_cast<int>(xoff);
    const auto yo = static_cast<int>(yoff);
    const auto zo = static_cast<int>(zoff);
    const std::array cidx0 = {cidx1[0] + xo, cidx1[1] + yo, cidx1[2] + zo};
    return {Segment{cidx0, p.old_location, pr0, true},
            Segment{cidx1, pr1, p.location, true}};
  }

  struct CurrentDeposition {
    using emdata_t = electromagnetics::EMData;
    using group_t = ParticleGroup;
    using CellData = ParticleGroup::CellData;
    using array_t = std::array<std::size_t, 3>;
    using EMFace = electromagnetics::EMFace;
    using EMSide = electromagnetics::EMSide;

    explicit CurrentDeposition(const emdata_t&) {}

    template<int D>
    static void updateJ(auto& J, const auto& as0, const auto& as1, const auto& bs0, const auto& bs1, const auto& cs0, const auto& cs1, const auto& qA, const std::array<std::size_t, 3>& idxs, const std::array<std::size_t, 3>& bounds) {
      static constexpr auto third = 1.0_fp / 3.0_fp;
      static constexpr auto sixth = 1.0_fp / 6.0_fp;

      const auto& [i, j, k] = idxs;
      const auto& [i1, j1, k1] = bounds;

      compute_t wm;
      auto wT = 0.0_fp;
      const std::array<compute_t, 2> ws = {as0[0] - as1[0], as1[2] - as0[2]};
      for (std::size_t ii = 0; ii < i1; ii++) {
        for (std::size_t jj = 0; jj < j1; jj++) {
          for (std::size_t kk = 0; kk < k1; kk++) {
            if constexpr      (D == 0) { wm = ws[ii]; }
            else if constexpr (D == 1) { wm = ws[jj]; }
            else                       { wm = ws[kk]; }
            wT = qA * (third * (bs0[jj] * cs0[kk] + bs1[jj] * cs1[kk]) + sixth * (bs1[jj] * cs0[kk] + bs0[jj] * cs1[kk]));
            // todo: need to be concerned about edges of grid here, where j - 1 < 0, j + 2 > Ny,... etc
#pragma omp atomic update
            J(i + ii, j + jj, k + kk) += wm * wT;
          }
        }
      }
    } // end updateJ()

    static void updateCell(const CellData& cell, const compute_t charge, emdata_t& emdata) {
      const auto& [ci, cj, ck] = cell.idxs;
      for (const auto& chunk : cell.chunks) {
        for (std::size_t pid = 0; pid < ParticleChunk::n_particles; pid++) {
          if (!chunk.active.test(pid)) { continue; }
          const auto& p = chunk[pid];

          for (const auto& [cids, p0, p1, active]: split_trajectory(p, {ci, cj, ck})) {
            if (!active) { continue; }
            const auto& [i, j, k] = cids;

            const auto xs0 = quad_shapes(p0[0]);
            const auto xs1 = quad_shapes(p1[0]);

            const auto ys0 = quad_shapes(p0[1]);
            const auto ys1 = quad_shapes(p1[1]);

            const auto zs0 = quad_shapes(p0[2]);
            const auto zs1 = quad_shapes(p1[2]);

            const auto cx = p.weight * charge / (Ayz * dt);
            const auto cy = p.weight * charge / (Axz * dt);
            const auto cz = p.weight * charge / (Axy * dt);

            updateJ<0>(emdata.Jx, xs0, xs1, ys0, ys1, zs0, zs1, cx, {i, j - 1, k - 1}, {2, 3, 3});
            updateJ<1>(emdata.Jy, ys0, ys1, xs0, xs1, zs0, zs1, cy, {i - 1, j, k - 1}, {3, 2, 3});
            updateJ<2>(emdata.Jz, zs0, zs1, xs0, xs1, ys0, ys1, cz, {i - 1, j - 1, k}, {3, 3, 2});
          } // end for(trajectory)
        } // end for(pid)
      } // end for(chunk)
    } // end update()

    static void update(const group_t& g, emdata_t& emdata) {
#pragma omp parallel for num_threads(nThreads) schedule(dynamic, chunkSize)
      for (std::size_t i = 0; i < g.cells.size(); i++) {
        if (g.cells[i].chunks.empty()) { continue; }
        updateCell(g.cells[i], g.charge, emdata);
      } // end for(i)
    }

    static void operator()(const group_t& g, emdata_t& emdata) {
      update(g, emdata);
    }
  }; // end struct CurrentDeposition
} // end namepsace tf::particles
#endif //CURRENT_DEPOSITION_HPP
