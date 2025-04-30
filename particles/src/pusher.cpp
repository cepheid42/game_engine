#include "pusher.hpp"

#include "program_params.hpp"
#include "em_params.hpp"
#include "constants.hpp"

#include <cassert>

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

  void BorisPush::update_velocity(Particle& p,
                                  const efield& Ex_c,
                                  const efield& Ey_c,
                                  const efield& Ez_c,
                                  const bfield& Bx_c,
                                  const bfield& By_c,
                                  const bfield& Bz_c)
  {
    const auto eps = EFieldAtParticle(p.location, Ex_c, Ey_c, Ez_c);
    const auto bet = BFieldAtParticle(p.location, Bx_c, By_c, Bz_c);

    const auto um = p.velocity.as_type<double>() + eps;
    const auto t = bet / std::sqrt(1.0 + um.length_squared() * constants::over_c_sqr<compute_t>);
    const auto s = 2.0 * t / (1.0 + t.length_squared());

    const auto v = eps + um + cross(um + cross(um, t), s);
    const auto gamma = std::sqrt(1.0 + v.length_squared() * constants::over_c_sqr<compute_t>);

    p.gamma = gamma;
    p.velocity = v.as_type<compute_t>();
  } // end BorisPush::update_velocity()

  bool BorisPush::update_position(Particle& p) {
    static constexpr vec3 delta_inv{dt / dx, dt / dy, dt / dz};

    const auto new_loc = p.location + delta_inv * p.velocity / p.gamma;
    const vec3 offsets = {std::floor(new_loc[0]), std::floor(new_loc[1]), std::floor(new_loc[2])};
    p.old_location = p.location - offsets;
    p.location = new_loc - offsets;
    return offsets[0] != 0 or offsets[2] != 0;
  } // end BorisPush::update_position()

  void BorisPush::advance_cell(cell_data& cell, const emdata_t& emdata, const compute_t qdt_over_2m) {
    const auto& [i, j, k] = cell.idxs;
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
      for (std::size_t p = 0; p < ParticleChunk::n_particles; p++) {
        if (!chunk.active.test(p)) { continue; }
        update_velocity(chunk[p], Ex_c, Ey_c, Ez_c, Bx_c, By_c, Bz_c);
      }
    }

    for (auto& chunk : cell.chunks) {
      for (std::size_t p = 0; p < ParticleChunk::n_particles; p++) {
        if (!chunk.active.test(p)) { continue; }
        chunk.moves.set(p, update_position(chunk[p]));
      }
    }
  } // end BorisPush::advance_cell()

  void BorisPush::advance(group_t& g, const emdata_t& emdata) {
#pragma omp parallel for num_threads(nThreads)
    for (std::size_t i = 0; i < Ncells; i++) {
      if (g.cells[i].chunks.empty()) { continue; }
      advance_cell(g.cells[i], emdata, g.qdt_over_2m);
    }
  } // end BorisPush::advance()

  void BorisPush::update_cells(group_t& group) {
    static constexpr std::size_t BC_DEPTH = PMLDepth - 5;
    for (auto& [chunks, idxs] : group.cells) {
      if (chunks.empty()) { continue; }
      const auto& [i, j, k] = idxs;
      for (auto& chunk : chunks) {
        for (std::size_t pid = 0; pid < ParticleChunk::n_particles; pid++) {
          if (!chunk.active.test(pid) or !chunk.moves.test(pid)) { continue; }
          auto& p = chunk[pid];

          const int x_offset = static_cast<int>(std::floor(p.old_location[0]));
          const int z_offset = static_cast<int>(std::floor(p.old_location[2]));
          const std::size_t inew = i - x_offset;
          const std::size_t knew = k - z_offset;

          // outflow boundary in X/Z
          if (inew < BC_DEPTH or inew > Nx - BC_DEPTH or knew < BC_DEPTH or knew > Nz - BC_DEPTH) {
            group.remove_particle(chunk, pid);
            continue;
          }

          buffer.emplace_back(p, get_cid(inew, j, knew));
          group.remove_particle(chunk, pid);
        } // end for(pid)
      } // end for(chunk)
    } // end for(cell)

    for (const auto& [p, cid]: buffer) {
      group.add_particle(p, cid);
    }
    buffer.clear();
  } // end BorisPush::update_cells()

  // void BorisPush::update_cells(group_t& group) {
  //   using buffer_t = std::vector<std::tuple<Particle, std::uint32_t, std::uint32_t, std::uint32_t>>;
  //   buffer_t buffer{};
  //   for (auto& [chunks, idxs] : group.cells) {
  //     const auto& [i, j, k] = idxs;
  //     for (auto& chunk : chunks) {
  //       for (std::size_t pid = 0; pid < ParticleChunk::n_particles; pid++) {
  //         if (!chunk.active.test(pid)) { continue; }
  //         auto& p = chunk[pid];
  //
  //         const int x_offset = static_cast<int>(std::floor(p.old_location[0]));
  //         // const int y_offset = static_cast<int>(std::floor(p.old_location[1]));
  //         const int z_offset = static_cast<int>(std::floor(p.old_location[2]));
  //
  //
  //         // Did not move out of the current cell
  //         if (x_offset == 0 and z_offset == 0) { continue; }
  //
  //         const std::size_t inew = i - x_offset;
  //         // const std::size_t jnew = j - y_offset;
  //         const std::size_t knew = k - z_offset;
  //
  //         // outflow boundary in X/Z
  //         if (inew < PMLDepth - 5 or inew > Nx - PMLDepth - 5 or knew < PMLDepth or knew > Nz - PMLDepth) {
  //           group.remove_particle(chunk, pid);
  //           continue;
  //         }
  //
  //         assert(get_cid(inew, j, knew) != get_cid(i, j, k));
  //
  //         // group.add_particle(p, inew, jnew, knew);
  //         buffer.emplace_back(p, inew, j, knew);
  //         group.remove_particle(chunk, pid);
  //       } // end for(pid)
  //     } // end for(chunk)
  //   } // end for(cell)
  //
  //   for (const auto& [p, i, j, k]: buffer) {
  //     group.add_particle(p, i, j, k);
  //     // group.add_particle(p, get_cid(i, j, k));
  //   }
  // } // end BorisPush::update_cells()
} // end namespace tf::particles