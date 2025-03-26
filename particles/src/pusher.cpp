#include "pusher.hpp"

#include "program_params.hpp"
#include "constants.hpp"

namespace tf::particles {
  vec3<compute_t> BFieldAtParticle(const compute_t fx, const compute_t fy, const compute_t fz,
                                   const std::array<compute_t, 2>& Hx_c,
                                   const std::array<compute_t, 2>& Hy_c,
                                   const std::array<compute_t, 2>& Hz_c)
  {
    return {(1.0f - fx) * Hx_c[0] + fx * Hx_c[1],
            (1.0f - fy) * Hy_c[0] + fy * Hy_c[1],
            (1.0f - fz) * Hz_c[0] + fz * Hz_c[1]};
  }

  vec3<compute_t> EFieldAtParticle(const compute_t fx, const compute_t fy, const compute_t fz,
                                   const std::array<compute_t, 4>& Ex_c,
                                   const std::array<compute_t, 4>& Ey_c,
                                   const std::array<compute_t, 4>& Ez_c)
  {
    const auto xy1 = (1.0f - fx) * (1.0f - fy);
    const auto xy2 = (1.0f - fx) * fy;
    const auto xy3 = fx * (1.0f - fy);
    const auto xy4 = fx * fy;

    const auto xz1 = (1.0f - fx) * (1.0f - fz);
    const auto xz2 = (1.0f - fx) * fz;
    const auto xz3 = fx * (1.0f - fz);
    const auto xz4 = fx * fz;

    const auto yz1 = (1.0f - fy) * (1.0f - fz);
    const auto yz2 = (1.0f - fy) * fz;
    const auto yz3 = fy * (1.0f - fz);
    const auto yz4 = fy * fz;

    return {yz1 * Ex_c[0] + yz2 * Ex_c[1] + yz3 * Ex_c[2] + yz4 * Ex_c[3],
            xz1 * Ey_c[0] + xz2 * Ey_c[1] + xz3 * Ey_c[2] + xz4 * Ey_c[3],
            xy1 + Ez_c[0] + xy2 + Ez_c[1] + xy3 + Ez_c[2] + xy4 * Ez_c[3]};
  }

  void BorisPush::update_particle(Particle& p,
                                  const ParticleGroup& g,
                                  const efield_comp& Ex_c,
                                  const efield_comp& Ey_c,
                                  const efield_comp& Ez_c,
                                  const bfield_comp& Bx_c,
                                  const bfield_comp& By_c,
                                  const bfield_comp& Bz_c)
  {
    const auto sx = p.location[0] / dx;
    const auto sy = p.location[1] / dy;
    const auto sz = p.location[2] / dz;

    const auto fx = sx - std::floor(sx);
    const auto fy = sy - std::floor(sy);
    const auto fz = sz - std::floor(sz);

    auto eps = g.qdt_over_2m * EFieldAtParticle(fx, fy, fz, Ex_c, Ey_c, Ez_c);
    const auto bet = g.qdt_over_2m * BFieldAtParticle(fx, fy, fz, Bx_c, By_c, Bz_c);

    const auto um = p.momentum * g.inv_mass + eps;
    const auto t = bet / static_cast<compute_t>(std::sqrt(1.0 + um.length_squared() * constants::over_c_sqr));
    const auto s = 2.0f * t / (1.0f + t.length_squared());

    eps += um + cross(um + cross(um, t), s);

    // Hopefully this give enough precision for calculating gamma
    const auto eps_double = math::SQR(static_cast<double>(eps[0]))
                          + math::SQR(static_cast<double>(eps[1]))
                          + math::SQR(static_cast<double>(eps[2]));

    p.momentum = g.mass * eps;
    p.gamma = std::sqrt(1.0 + eps_double * constants::over_c_sqr);
    // todo: updating location here, since it seems reasonable
    p.old_location = p.location;
    p.location += dt * g.inv_mass * p.momentum / static_cast<compute_t>(p.gamma);
  }

  void BorisPush::update_cell(p_vec& particles, const std::array<std::size_t, 3>& coords, const group_t& g, const emdata_t& emdata) {
    const auto& [i, j, k] = coords;

    const std::array Ex_c = {emdata.Ex(i, j, k), emdata.Ex(i, j + 1, k), emdata.Ex(i, j, k + 1), emdata.Ex(i, j + 1, k + 1)};
    const std::array Ey_c = {emdata.Ey(i, j, k), emdata.Ey(i + 1, j, k), emdata.Ey(i, j, k + 1), emdata.Ey(i + 1, j, k + 1)};
    const std::array Ez_c = {emdata.Ez(i, j, k), emdata.Ez(i, j + 1, k), emdata.Ez(i + 1, j, k), emdata.Ez(i + 1, j + 1, k)};

    // todo: B-fields are not aligned in time with E-fields at this point
    const std::array Bx_c = {emdata.Bx(i, j, k), emdata.Bx(i + 1, j, k)};
    const std::array By_c = {emdata.By(i, j, k), emdata.By(i, j + 1, k)};
    const std::array Bz_c = {emdata.Bz(i, j, k), emdata.Bz(i, j, k + 1)};

    for (std::size_t cid = 0; cid < particles.size(); cid++) {
      update_particle(particles[cid], g, Ex_c, Ey_c, Ez_c, Bx_c, By_c, Bz_c);
    }
  }

  void BorisPush::visit(const p_tree& node, const group_t& g, const emdata_t& emdata) {
    for (std::size_t i = 0; i < 8; i++) {
      if (node.is_leaf) {
        if (node.active.test(i)) {
          update_cell(*node.cells[i], node.cell_coords, g, emdata);
        }
      } else {
        visit(node.children[i], g, emdata);
      }
    }
  }
}