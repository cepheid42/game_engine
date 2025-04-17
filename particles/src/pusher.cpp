#include "pusher.hpp"

#include "program_params.hpp"
#include "constants.hpp"

#define UTL_PROFILER_DISABLE
#include "profiler.hpp"

#include <cassert>

#include "em_params.hpp"

namespace tf::particles {
  vec3<double> BFieldAtParticle(const vec3<compute_t>& loc,
                                const std::array<compute_t, 2>& Bx_c,
                                const std::array<compute_t, 2>& By_c,
                                const std::array<compute_t, 2>& Bz_c)
  {
    return {(1.0 - loc[0]) * Bx_c[0] + loc[0] * Bx_c[1],
            (1.0 - loc[1]) * By_c[0] + loc[1] * By_c[1],
            (1.0 - loc[2]) * Bz_c[0] + loc[2] * Bz_c[1]};
  } // end BFieldAtParticle()

  vec3<double> EFieldAtParticle(const vec3<compute_t>& loc,
                                const std::array<compute_t, 4>& Ex_c,
                                const std::array<compute_t, 4>& Ey_c,
                                const std::array<compute_t, 4>& Ez_c)
  {
    const auto xy1 = (1.0 - loc[0]) * (1.0 - loc[1]);
    const auto xy2 = (1.0 - loc[0]) * loc[1];
    const auto xy3 = loc[0] * (1.0 - loc[1]);
    const auto xy4 = loc[0] * loc[1];

    const auto xz1 = (1.0 - loc[0]) * (1.0 - loc[2]);
    const auto xz2 = (1.0 - loc[0]) * loc[2];
    const auto xz3 = loc[0] * (1.0 - loc[2]);
    const auto xz4 = loc[0] * loc[2];

    const auto yz1 = (1.0 - loc[1]) * (1.0 - loc[2]);
    const auto yz2 = (1.0 - loc[1]) * loc[2];
    const auto yz3 = loc[1] * (1.0 - loc[2]);
    const auto yz4 = loc[1] * loc[2];

    return {yz1 * Ex_c[0] + yz2 * Ex_c[1] + yz3 * Ex_c[2] + yz4 * Ex_c[3],
            xz1 * Ey_c[0] + xz2 * Ey_c[1] + xz3 * Ey_c[2] + xz4 * Ey_c[3],
            xy1 * Ez_c[0] + xy2 * Ez_c[1] + xy3 * Ez_c[2] + xy4 * Ez_c[3]};
  } // end EFieldAtParticle()

  void BorisPush::update_particle(Particle& p,
                                  const efield& Ex_c,
                                  const efield& Ey_c,
                                  const efield& Ez_c,
                                  const bfield& Bx_c,
                                  const bfield& By_c,
                                  const bfield& Bz_c)
  {
    // Inverse square root intrinsics
    // https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=rsqrt&expand=4804&ig_expand=5653,5653
    constexpr vec3 delta_inv{1.0 / dx, 1.0 / dy, 1.0 / dz};

    const auto eps = EFieldAtParticle(p.location, Ex_c, Ey_c, Ez_c);
    const auto bet = BFieldAtParticle(p.location, Bx_c, By_c, Bz_c);

    const auto um = p.velocity.as_type<double>() + eps;
    const auto t = bet / std::sqrt(1.0 + um.as_type<double>().length_squared() * constants::over_c_sqr);
    const auto s = 2.0 * t / (1.0 + t.length_squared());

    // const auto v = eps + um + cross(um + cross(um, t), s);
    auto v = eps + um + cross_simd_dub(um + cross_simd_dub(um, t), s);
    v[1] = 0.0_fp;
    const auto gamma = std::sqrt(1.0 + v.length_squared() * constants::over_c_sqr);

    // todo: updating location here, since it seems reasonable
    const auto new_loc = p.location + ((static_cast<double>(dt) / gamma) * (v * delta_inv)).as_type<compute_t>();
    const vec3 offsets = {std::floor(new_loc[0]), std::floor(new_loc[1]), std::floor(new_loc[2])};

    p.velocity = v.as_type<compute_t>();
    p.old_location = p.location - offsets;
    p.location = new_loc - offsets;
  } // end BorisPush::update_particle()

  void BorisPush::advance_cell(cell_data& cell, const emdata_t& emdata, const compute_t qdt_over_2m, const std::size_t i, const std::size_t j, const std::size_t k) {
    efield Ex_c = {emdata.getEx(i, j, k), emdata.getEx(i, j + 1, k), emdata.getEx(i, j, k + 1), emdata.getEx(i, j + 1, k + 1)};
    efield Ey_c = {emdata.getEy(i, j, k), emdata.getEy(i + 1, j, k), emdata.getEy(i, j, k + 1), emdata.getEy(i + 1, j, k + 1)};
    efield Ez_c = {emdata.getEz(i, j, k), emdata.getEz(i, j + 1, k), emdata.getEz(i + 1, j, k), emdata.getEz(i + 1, j + 1, k)};
    for (std::size_t q = 0; q < 4; q++) {
      Ex_c[q] *= qdt_over_2m;
      Ey_c[q] *= qdt_over_2m;
      Ez_c[q] *= qdt_over_2m;
    }

    // todo: B-fields are not aligned in time with E-fields at this point
    bfield Bx_c = {emdata.getBx(i, j, k), emdata.getBx(i + 1, j, k)};
    bfield By_c = {emdata.getBy(i, j, k), emdata.getBy(i, j + 1, k)};
    bfield Bz_c = {emdata.getBz(i, j, k), emdata.getBz(i, j, k + 1)};
    for (std::size_t q = 0; q < 2; q++) {
      Bx_c[q] *= qdt_over_2m;
      By_c[q] *= qdt_over_2m;
      Bz_c[q] *= qdt_over_2m;
    }

    for (auto& chunk : cell.chunks) {
#pragma omp simd
      for (std::size_t p = 0; p < ParticleChunk::n_particles; p++) {
        update_particle(chunk[p], Ex_c, Ey_c, Ez_c, Bx_c, By_c, Bz_c);
      }
    }
  } // end BorisPush::advance_particles()

  void BorisPush::advance(group_t& g, const emdata_t& emdata) {
#pragma omp parallel for collapse(3) num_threads(nThreads)
    for (std::size_t i = 0; i < Ncx; i++) {
      for (std::size_t j = 0; j < Ncy; j++) {
        for (std::size_t k = 0; k < Ncz; k++) {
          const auto index = get_cid(i, j, k);
          advance_cell(g.cells[index], emdata, g.qdt_over_2m, i, j, k);
        } // end for(k)
      } // end for(j)
    } // end for(i)
  } // end BorisPush::visit()

  void BorisPush::update_cells(group_t& group) {
    // constexpr std::size_t alpha_x = (Ncx - 1) - (nHalo + 1);
    // constexpr std::size_t alpha_y = (Ncy - 1) - (nHalo + 1);
    // constexpr std::size_t alpha_z = (Ncz - 1) - (nHalo + 1);

    using buffer_t = std::vector<std::tuple<Particle, std::uint32_t, std::uint32_t, std::uint32_t>>;
    buffer_t buffer{};
    for (auto& [chunks, idxs] : group.cells) {
      const auto& [i, j, k] = idxs;
      for (auto& chunk : chunks) {
        for (std::size_t pid = 0; pid < ParticleChunk::n_particles; pid++) {
          if (!chunk.active.test(pid)) { continue; }
          auto& p = chunk[pid];

          const int x_offset = static_cast<int>(std::floor(p.old_location[0]));
          const int y_offset = static_cast<int>(std::floor(p.old_location[1]));
          const int z_offset = static_cast<int>(std::floor(p.old_location[2]));

          assert(y_offset == 0);

          // Did not move out of the current cell
          if (x_offset == 0 and y_offset == 0 and z_offset == 0) { continue; }

          const std::size_t inew = i - x_offset;
          const std::size_t jnew = j - y_offset;
          const std::size_t knew = k - z_offset;

          // outflow boundary in X/Z
          if (inew < PMLDepth - 5 or inew > Nx - PMLDepth - 5 or knew < PMLDepth or knew > Nz - PMLDepth) {
            group.remove_particle(chunk, pid);
            continue;
          }

          // // Periodic in X
          // std::size_t inew = i - x_offset;
          // if (inew < nHalo) {
          //   inew += alpha_x;
          // } else if (inew > (Ncx - nHalo - 1)) {
          //   inew -= alpha_x;
          // }
          //
          // // Periodic in Y
          // std::size_t jnew = j - y_offset;
          // if (jnew < nHalo) {
          //   jnew += alpha_y;
          // } else if (jnew > (Ncy - nHalo - 1)) {
          //   jnew -= alpha_y;
          // }
          //
          // // Periodic in Z
          // std::size_t knew = k - z_offset;
          // if (knew < nHalo) {
          //   knew += alpha_z;
          // } else if (knew > (Ncz - nHalo - 1)) {
          //   knew -= alpha_z;
          // }

          assert(get_cid(inew, jnew, knew) != get_cid(i, j, k));

          // group.add_particle(p, inew, jnew, knew);
          buffer.emplace_back(p, inew, jnew, knew);
          group.remove_particle(chunk, pid);
        } // end for(pid)
      } // end for(chunk)
    } // end for(cell)

    for (const auto& [p, i, j, k]: buffer) {
      group.add_particle(p, i, j, k);
    }
  } // end BorisPush::update_cells()
} // end namespace tf::particles