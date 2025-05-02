#include "pusher.hpp"

#include "program_params.hpp"
#include "em_params.hpp"
#include "constants.hpp"

#include <cassert>

namespace tf::particles {
  std::array<double, 2> shapesCIC(const auto x) { return {1.0 - x, x}; }

  std::array<double, 3> shapesQuad(const auto x) {
    return {0.5 * math::SQR(0.5 - x),
            0.75 - math::SQR(x),
            0.5 * math::SQR(0.5 + x)};
  }

  template<std::size_t nx, std::size_t ny, std::size_t nz>
  std::array<double, nx*ny*nz> computeTotalShape(const auto& xs, const auto& ys, const auto& zs) {
    std::array<double, nx*ny*nz> shapes{};
    for (std::size_t i = 0; i < nx; i++) {
      for (std::size_t j = 0; j < ny; j++) {
        for (std::size_t k = 0; k < nz; k++) {
          const auto idx = k + (nz * j) + (ny * nz * i);
          shapes[idx] = xs[i] * ys[j] * zs[k];
        }
      }
    }
    return shapes;
  }

  template<std::size_t F, std::size_t nx, std::size_t ny, std::size_t nz>
  auto getField(const auto& emdata, const auto& idxs, const std::array<int, 6>& offsets) {
    std::array<double, nx*ny*nz> fields{};
    int idx = 0;
    for (int i = offsets[0]; i <= offsets[1]; i++) {
      for (int j = offsets[2]; j < offsets[3]; j++) { // should only run once
        for (int k = offsets[4]; k <= offsets[5]; k++) {
          fields[idx] = emdata.template get<F>(idxs[0] + i, idxs[1] + j, idxs[2] + k);
          idx++;
        }
      }
    }
    return fields;
  }

  template<std::size_t N>
  auto calculateShape(const auto& field, const auto& shape) -> double {
    double result = 0.0;
    for (std::size_t i = 0; i < N; i++) {
      result += field[i] * shape[i];
    }
    return result;
  }

  std::array<double, 6> FieldAtParticle(Particle& p,
                                         const auto& Ex_c,
                                         const auto& Ey_c,
                                         const auto& Ez_c,
                                         const auto& Bx_c,
                                         const auto& By_c,
                                         const auto& Bz_c)
  {
    const auto xr_shapes = shapesCIC(p.location[0]);
    const auto yr_shapes = shapesCIC(p.location[1]);
    const auto zr_shapes = shapesCIC(p.location[2]);
    const auto x_shapes = shapesQuad(p.location[0]);
    const auto y_shapes = shapesQuad(p.location[1]);
    const auto z_shapes = shapesQuad(p.location[2]);

    const auto ex_shapes = computeTotalShape<2, 1, 3>(xr_shapes, y_shapes, z_shapes);
    const auto ey_shapes = computeTotalShape<3, 1, 3>(x_shapes, yr_shapes, z_shapes);
    const auto ez_shapes = computeTotalShape<3, 1, 2>(x_shapes, y_shapes, zr_shapes);
    const auto bx_shapes = computeTotalShape<3, 1, 2>(x_shapes, yr_shapes, zr_shapes);
    const auto by_shapes = computeTotalShape<2, 1, 2>(xr_shapes, y_shapes, zr_shapes);
    const auto bz_shapes = computeTotalShape<2, 1, 3>(xr_shapes, yr_shapes, z_shapes);

    return {calculateShape<6>(Ex_c, ex_shapes),
            calculateShape<9>(Ey_c, ey_shapes),
            calculateShape<6>(Ez_c, ez_shapes),
            calculateShape<6>(Bx_c, bx_shapes),
            calculateShape<4>(By_c, by_shapes),
            calculateShape<6>(Bz_c, bz_shapes)};
  }



  void BorisPush::update_velocity(Particle& p,
                                  const auto& Ex_c,
                                  const auto& Ey_c,
                                  const auto& Ez_c,
                                  const auto& Bx_c,
                                  const auto& By_c,
                                  const auto& Bz_c,
                                  const auto qdt)
  {
    const auto emf = FieldAtParticle(p, Ex_c, Ey_c, Ez_c, Bx_c, By_c, Bz_c);

    const vec3<double> eps{qdt * emf[0], qdt * emf[1], qdt * emf[2]};
    const vec3<double> bet{qdt * emf[3], qdt * emf[4], qdt * emf[5]};

    const auto um = p.velocity.as_type<double>() + eps;
    const auto t = bet / std::sqrt(1.0 + um.length_squared() * constants::over_c_sqr<double>);
    const auto s = 2.0 * t / (1.0 + t.length_squared());

    const auto v = eps + um + cross(um + cross(um, t), s);
    const auto gamma = std::sqrt(1.0 + v.length_squared() * constants::over_c_sqr<double>);

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
    const auto Ex_c = getField<0,2,1,3>(emdata, cell.idxs, { 0, 1, 0, 1, -1, 1});
    const auto Ey_c = getField<1,3,1,3>(emdata, cell.idxs, {-1, 1, 0, 1, -1, 1});
    const auto Ez_c = getField<2,3,1,2>(emdata, cell.idxs, {-1, 1, 0, 1,  0, 1});

    const auto Bx_c = getField<3,3,1,2>(emdata, cell.idxs, {-1, 1, 0, 1,  0, 1});
    const auto By_c = getField<4,2,1,2>(emdata, cell.idxs, { 0, 1, 0, 1,  0, 1});
    const auto Bz_c = getField<5,2,1,3>(emdata, cell.idxs, { 0, 1, 0, 1, -1, 1});

    for (auto& chunk : cell.chunks) {
      for (std::size_t p = 0; p < ParticleChunk::n_particles; p++) {
        if (!chunk.active.test(p)) { continue; }
        update_velocity(chunk[p], Ex_c, Ey_c, Ez_c, Bx_c, By_c, Bz_c, qdt_over_2m);
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
#pragma omp parallel for collapse(3) num_threads(nThreads) schedule(dynamic, chunkSize)
    for (std::size_t i = BC_DEPTH; i < Ncx - BC_DEPTH; i++) {
      for (std::size_t j = 0; j < Ncy; j++) {
        for (std::size_t k = BC_DEPTH; k < Ncz - BC_DEPTH; k++) {
          const auto idx = get_cid(i, j, k);
          if (g.cells[idx].chunks.empty()) { continue; }
          advance_cell(g.cells[idx], emdata, g.qdt_over_2m);
        }
      }
    }
  } // end BorisPush::advance()

  void BorisPush::update_cells(group_t& group) {
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
} // end namespace tf::particles