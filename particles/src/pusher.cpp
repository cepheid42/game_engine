#include "pusher.hpp"

#include "program_params.hpp"
#include "constants.hpp"

namespace tf::particles {
  vec3<compute_t> BFieldAtParticle(const vec3<compute_t>& loc,
                                   const std::array<compute_t, 2>& Hx_c,
                                   const std::array<compute_t, 2>& Hy_c,
                                   const std::array<compute_t, 2>& Hz_c)
  {
    return {(1.0f - loc[0]) * Hx_c[0] + loc[0] * Hx_c[1],
            (1.0f - loc[1]) * Hy_c[0] + loc[1] * Hy_c[1],
            (1.0f - loc[2]) * Hz_c[0] + loc[2] * Hz_c[1]};
  }

  vec3<compute_t> EFieldAtParticle(const vec3<compute_t>& loc,
                                   const std::array<compute_t, 4>& Ex_c,
                                   const std::array<compute_t, 4>& Ey_c,
                                   const std::array<compute_t, 4>& Ez_c)
  {
    const auto xy1 = (1.0f - loc[0]) * (1.0f - loc[1]);
    const auto xy2 = (1.0f - loc[0]) * loc[1];
    const auto xy3 = loc[0] * (1.0f - loc[1]);
    const auto xy4 = loc[0] * loc[1];

    const auto xz1 = (1.0f - loc[0]) * (1.0f - loc[2]);
    const auto xz2 = (1.0f - loc[0]) * loc[2];
    const auto xz3 = loc[0] * (1.0f - loc[2]);
    const auto xz4 = loc[0] * loc[2];

    const auto yz1 = (1.0f - loc[1]) * (1.0f - loc[2]);
    const auto yz2 = (1.0f - loc[1]) * loc[2];
    const auto yz3 = loc[1] * (1.0f - loc[2]);
    const auto yz4 = loc[1] * loc[2];

    return {yz1 * Ex_c[0] + yz2 * Ex_c[1] + yz3 * Ex_c[2] + yz4 * Ex_c[3],
            xz1 * Ey_c[0] + xz2 * Ey_c[1] + xz3 * Ey_c[2] + xz4 * Ey_c[3],
            xy1 + Ez_c[0] + xy2 + Ez_c[1] + xy3 + Ez_c[2] + xy4 * Ez_c[3]};
  }

  void BorisPush::update_particle(Particle& p,
                                  const compute_t qdt_over_2m,
                                  const efield& Ex_c,
                                  const efield& Ey_c,
                                  const efield& Ez_c,
                                  const bfield& Bx_c,
                                  const bfield& By_c,
                                  const bfield& Bz_c)
  {
    auto eps = qdt_over_2m * EFieldAtParticle(p.location, Ex_c, Ey_c, Ez_c);
    const auto bet = qdt_over_2m * BFieldAtParticle(p.location, Bx_c, By_c, Bz_c);

    const auto um = p.velocity + eps;
    const auto t = bet / static_cast<compute_t>(std::sqrt(1.0 + um.length_squared() * constants::over_c_sqr));
    const auto s = 2.0f * t / (1.0f + t.length_squared());

    eps += um + cross(um + cross(um, t), s);

    p.velocity = eps.to_float();
    p.gamma = std::sqrt(1.0 + eps.length_squared() * constants::over_c_sqr);
    
    // todo: updating location here, since it seems reasonable
    const auto new_loc = p.location + (dt * p.velocity / static_cast<compute_t>(p.gamma));
    const vec3<compute_t> offsets = {std::floor(new_loc[0]), std::floor(new_loc[1]), std::floor(new_loc[2])};

    p.old_location = p.location - offsets;
    p.location = new_loc - offsets;
  }

  void BorisPush::advance_particles(ParticleCell& c, const group_t& g, const emdata_t& emdata) {
    auto& particles = c.particles;
    const auto& [i, j, k] = morton_decode(c.cid);

    const efield Ex_c = {emdata.Ex(i, j, k), emdata.Ex(i, j + 1, k), emdata.Ex(i, j, k + 1), emdata.Ex(i, j + 1, k + 1)};
    const efield Ey_c = {emdata.Ey(i, j, k), emdata.Ey(i + 1, j, k), emdata.Ey(i, j, k + 1), emdata.Ey(i + 1, j, k + 1)};
    const efield Ez_c = {emdata.Ez(i, j, k), emdata.Ez(i, j + 1, k), emdata.Ez(i + 1, j, k), emdata.Ez(i + 1, j + 1, k)};

    // todo: B-fields are not aligned in time with E-fields at this point
    const bfield Bx_c = {emdata.Bx(i, j, k), emdata.Bx(i + 1, j, k)};
    const bfield By_c = {emdata.By(i, j, k), emdata.By(i, j + 1, k)};
    const bfield Bz_c = {emdata.Bz(i, j, k), emdata.Bz(i, j, k + 1)};

    for (std::size_t cid = 0; cid < particles.size(); cid++) {
      update_particle(particles[cid], g.qdt_over_2m, Ex_c, Ey_c, Ez_c, Bx_c, By_c, Bz_c);
    }
  } // end BorisPush::advance_particles()

  void BorisPush::push_particles(const p_tree& node, const group_t& g, const emdata_t& emdata) {
    for (std::size_t i = 0; i < 8; i++) {
      if (node.is_leaf) {
        if (node.active.test(i)) {
          advance_particles(*node.cells[i], g, emdata);
        }
      } else {
        push_particles(node.children[i], g, emdata);
      }
    }
  } // end BorisPush::visit()

  void BorisPush::update_cells(const p_tree& node, group_t& g) {
    // todo: Got to take care of over/underflow of the edges cells
    for (std::size_t c = 0; c < 8; c++) {
      if (node.is_leaf) {
        if (node.active.test(c)) {
          const auto& [i, j, k] = morton_decode(node.cells[c]->cid);
          // copy into new cells
          for (auto& p : node.cells[c]->particles) {
            const int x_offset = static_cast<int>(std::floor(p.old_location[0]));
            const int y_offset = static_cast<int>(std::floor(p.old_location[1]));
            const int z_offset = static_cast<int>(std::floor(p.old_location[2]));

            if (x_offset == 0 or y_offset == 0 or z_offset == 0) { continue; }

            const auto new_cid = morton_encode(i - x_offset, j - y_offset, k - z_offset);
            if (node.cells[c]->cid != new_cid) {
              g.cells[new_cid].particles.push_back(p);
              p.active = false;
            }
          }
          // erase all inactive particles from current cell
          std::erase_if(node.cells[c]->particles, [&](const auto& p) { return p.active; });
        }
      } else {
        update_cells(node.children[c], g);
      }
    }
  }

} // end namespace tf::particles