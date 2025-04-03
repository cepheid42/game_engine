#ifndef CURRENT_DEPOSITION_HPP
#define CURRENT_DEPOSITION_HPP

#include "math_utils.hpp"

#include "particles.hpp"
#include "em_data.hpp"
#include "bc_data.hpp"
#include "bc_functors.hpp"

#include <array>
#include <cmath>

namespace tf::particles {
  inline std::array<compute_t, 3> quad_shapes(const compute_t v) {
    return {
      0.5f * math::SQR(0.5f - v),
      0.75f - math::SQR(v),
      0.5f * math::SQR(0.5f + v)
    };
  }

  struct Segment {
    std::array<std::size_t, 3> cids = {0, 0, 0};
    vec3<compute_t> p0 = {0.0_fp, 0.0_fp, 0.0_fp};
    vec3<compute_t> p1 = {0.0_fp, 0.0_fp, 0.0_fp};
    bool active = false;
  };

  inline std::array<Segment, 2> split_trajectory(const Particle& p, const std::array<std::size_t, 3>& cidx1) {
    const auto x_offset = std::floor(p.old_location[0]);
    const auto y_offset = std::floor(p.old_location[1]);
    const auto z_offset = std::floor(p.old_location[2]);

    if (x_offset == 0.0 and y_offset == 0.0 and z_offset == 0.0) {
      // Particle did not leave current cell
      return {Segment{cidx1, p.old_location, p.location, true},
              Segment{cidx1, p.old_location, p.location, false}};
    }

    const auto xr0 = x_offset == 0.0_fp ? p.location[0] : std::max(0.0_fp, -x_offset);
    const auto yr0 = y_offset == 0.0_fp ? p.location[1] : std::max(0.0_fp, -y_offset);
    const auto zr0 = z_offset == 0.0_fp ? p.location[2] : std::max(0.0_fp, -z_offset);

    const vec3<compute_t> pr0{xr0, yr0, zr0};
    const vec3<compute_t> pr1{xr0 + x_offset, yr0 + y_offset, zr0 + z_offset};

    // todo: pray to the gods that none of the cell id's are zero...
    const auto xo = static_cast<int>(x_offset);
    const auto yo = static_cast<int>(y_offset);
    const auto zo = static_cast<int>(z_offset);
    const std::array<std::size_t, 3> cidx0 = {cidx1[0] + xo, cidx1[1] + yo, cidx1[2] + zo};
    return {Segment{cidx0, p.old_location, pr0, true},
            Segment{cidx1, pr1, p.location, true}};
  }

  template<electromagnetics::EMFace F>
  struct PeriodicBC {
    using PeriodicData = electromagnetics::PeriodicData<F, electromagnetics::EMSide::Lo>;
    using PeriodicFunctor = electromagnetics::PeriodicFunctor<F, true>;
    using BCIntegrator = electromagnetics::BCIntegrator<PeriodicFunctor>;

    static constexpr Array3D<void> empty{};

    BCIntegrator updater_func{};
    PeriodicData J0;
    PeriodicData J1;

    explicit PeriodicBC(const auto& j0, const auto& j1)
    : J0(j0), J1(j1)
    {}

    void operator()(auto& j0, auto& j1) {
      updater_func(j0, empty, empty, J0);
      updater_func(j1, empty, empty, J1);
    }
  };


  struct CurrentDeposition {
    using emdata_t = electromagnetics::EMData;
    using group_t = ParticleGroup;;
    using p_tree = Octree<ParticleCell>;
    using array_t = std::array<std::size_t, 3>;
    using EMFace = electromagnetics::EMFace;
    using EMSide = electromagnetics::EMSide;

    PeriodicBC<EMFace::X> x_bc;
    PeriodicBC<EMFace::Y> y_bc;
    PeriodicBC<EMFace::Z> z_bc;

    explicit CurrentDeposition(const emdata_t& emdata)
    : x_bc(emdata.Jy, emdata.Jz),
      y_bc(emdata.Jx, emdata.Jz),
      z_bc(emdata.Jx, emdata.Jy)
    {}

    template<int D>
    static void updateJ(auto& J, const auto& as0, const auto& as1, const auto& bs0, const auto& bs1, const auto& cs0, const auto& cs1, const auto& qA, const std::array<std::size_t, 3>& idxs, const std::array<std::size_t, 3>& bounds) {
      static constexpr auto third = 1.0f / 3.0f;
      static constexpr auto sixth = 1.0f / 6.0f;

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
            J(i + ii, j + jj, k + kk) += wm * wT;
          }
        }
      }
    } // end updateJ()

    static void updateCell(const std::vector<Particle>& particles, const std::size_t cid, const compute_t charge, emdata_t& emdata) {
      const auto& [i, j, k] = morton_decode(cid);
      for (std::size_t pid = 0; pid < particles.size(); pid++) {
        const auto& p = particles[pid];

        for (const auto& segment: split_trajectory(p, {i, j, k})) {
          if (!segment.active) { continue; }
          const auto& [i, j, k] = segment.cids;

          const auto xs0 = quad_shapes(segment.p0[0]);
          const auto xs1 = quad_shapes(segment.p1[0]);

          const auto ys0 = quad_shapes(segment.p0[1]);
          const auto ys1 = quad_shapes(segment.p1[1]);

          const auto zs0 = quad_shapes(segment.p0[2]);
          const auto zs1 = quad_shapes(segment.p1[2]);

          const auto cx = p.weight * charge / (Ayz * dt);
          const auto cy = p.weight * charge / (Axz * dt);
          const auto cz = p.weight * charge / (Axy * dt);

          updateJ<0>(emdata.Jx, xs0, xs1, ys0, ys1, zs0, zs1, cx, {i, j - 1, k - 1}, {2, 3, 3});
          updateJ<1>(emdata.Jy, ys0, ys1, xs0, xs1, zs0, zs1, cy, {i - 1, j, k - 1}, {3, 2, 3});
          updateJ<2>(emdata.Jz, zs0, zs1, xs0, xs1, ys0, ys1, cz, {i - 1, j - 1, k}, {3, 3, 2});
        }
      }
    } // end update()

    static void update(const p_tree& node, const group_t& g, emdata_t& emdata) {
      for (std::size_t i = 0; i < 8; i++) {
        if (!node.active.test(i)) { continue; }
        if (node.is_leaf) {
          updateCell(node.cells[i]->particles, node.cells[i]->cid, g.charge, emdata);
        } else {
          update(node.children[i], g, emdata);
        }
      }
    }

    void operator()(const group_t& g, emdata_t& emdata) {
      update(g.tree, g, emdata);
      x_bc(emdata.Jy, emdata.Jz);
      y_bc(emdata.Jx, emdata.Jz);
      z_bc(emdata.Jx, emdata.Jy);
    }
  }; // end struct CurrentDeposition
} // end namepsace tf::particles
#endif //CURRENT_DEPOSITION_HPP
