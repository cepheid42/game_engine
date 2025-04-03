#include "pusher.hpp"

#include "program_params.hpp"
#include "constants.hpp"

#include <cassert>

namespace tf::particles {
  vec3<compute_t> BFieldAtParticle(const vec3<compute_t>& loc,
                                   const std::array<compute_t, 2>& Bx_c,
                                   const std::array<compute_t, 2>& By_c,
                                   const std::array<compute_t, 2>& Bz_c)
  {
    return {(1.0_fp - loc[0]) * Bx_c[0] + loc[0] * Bx_c[1],
            (1.0_fp - loc[1]) * By_c[0] + loc[1] * By_c[1],
            (1.0_fp - loc[2]) * Bz_c[0] + loc[2] * Bz_c[1]};
  } // end BFieldAtParticle()

  vec3<compute_t> EFieldAtParticle(const vec3<compute_t>& loc,
                                   const std::array<compute_t, 4>& Ex_c,
                                   const std::array<compute_t, 4>& Ey_c,
                                   const std::array<compute_t, 4>& Ez_c)
  {
    const auto xy1 = (1.0_fp - loc[0]) * (1.0_fp - loc[1]);
    const auto xy2 = (1.0_fp - loc[0]) * loc[1];
    const auto xy3 = loc[0] * (1.0_fp - loc[1]);
    const auto xy4 = loc[0] * loc[1];

    const auto xz1 = (1.0_fp - loc[0]) * (1.0_fp - loc[2]);
    const auto xz2 = (1.0_fp - loc[0]) * loc[2];
    const auto xz3 = loc[0] * (1.0_fp - loc[2]);
    const auto xz4 = loc[0] * loc[2];

    const auto yz1 = (1.0_fp - loc[1]) * (1.0_fp - loc[2]);
    const auto yz2 = (1.0_fp - loc[1]) * loc[2];
    const auto yz3 = loc[1] * (1.0_fp - loc[2]);
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

    constexpr vec3 delta_inv{1.0_fp / dx, 1.0_fp / dy, 1.0_fp / dz};
    
    const auto eps = EFieldAtParticle(p.location, Ex_c, Ey_c, Ez_c);
    const auto bet = BFieldAtParticle(p.location, Bx_c, By_c, Bz_c);
    
    const auto um = p.velocity + eps;
    const auto t = bet / static_cast<compute_t>(std::sqrt(1.0 + um.as_type<double>().length_squared() * constants::over_c_sqr));
    const auto s = 2.0_fp * t / (1.0_fp + t.length_squared());

    p.velocity = eps + um + cross(um + cross(um, t), s);
    const auto gamma = std::sqrt(1.0 + p.velocity.as_type<double>().length_squared() * constants::over_c_sqr);

    // todo: updating location here, since it seems reasonable
    const auto new_loc = p.location + static_cast<compute_t>(static_cast<double>(dt) / gamma) * (p.velocity * delta_inv);
    const vec3<compute_t> offsets = {std::floor(new_loc[0]), std::floor(new_loc[1]), std::floor(new_loc[2])};

    p.old_location = p.location - offsets;
    p.location = new_loc - offsets;
  } // end BorisPush::update_particle()

  void BorisPush::advance_particles(ParticleCell& c, const group_t& g, const emdata_t& emdata) {
    auto& particles = c.particles;
    const auto& [i, j, k] = morton_decode(c.cid);

    efield Ex_c = {emdata.getEx(i, j, k), emdata.getEx(i, j + 1, k), emdata.getEx(i, j, k + 1), emdata.getEx(i, j + 1, k + 1)};
    efield Ey_c = {emdata.getEy(i, j, k), emdata.getEy(i + 1, j, k), emdata.getEy(i, j, k + 1), emdata.getEy(i + 1, j, k + 1)};
    efield Ez_c = {emdata.getEz(i, j, k), emdata.getEz(i, j + 1, k), emdata.getEz(i + 1, j, k), emdata.getEz(i + 1, j + 1, k)};
    for (std::size_t q = 0; q < 4; q++) {
      Ex_c[q] *= g.qdt_over_2m;
      Ey_c[q] *= g.qdt_over_2m;
      Ez_c[q] *= g.qdt_over_2m;
    }

    // todo: B-fields are not aligned in time with E-fields at this point
    bfield Bx_c = {emdata.getBx(i, j, k), emdata.getBx(i + 1, j, k)};
    bfield By_c = {emdata.getBy(i, j, k), emdata.getBy(i, j + 1, k)};
    bfield Bz_c = {emdata.getBz(i, j, k), emdata.getBz(i, j, k + 1)};
    for (std::size_t q = 0; q < 2; q++) {
      Bx_c[q] *= g.qdt_over_2m;
      By_c[q] *= g.qdt_over_2m;
      Bz_c[q] *= g.qdt_over_2m;
    }

    for (auto& particle : particles) {
      update_particle(particle, Ex_c, Ey_c, Ez_c, Bx_c, By_c, Bz_c);
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

            if (x_offset == 0 and y_offset == 0 and z_offset == 0) { continue; }

            auto i_new = i - x_offset;
            auto j_new = j - y_offset;
            auto k_new = k - z_offset;

            if (i_new == nHalo or i_new == Nx - nHalo) {
              i_new = Nx - i_new;
            }
            if (j_new == nHalo or j_new == Nx - nHalo) {
              j_new = Ny - j_new;
            }
            if (k_new == nHalo or k_new == Nx - nHalo) {
              k_new = Nz - k_new;
            }

            const auto new_cid = morton_encode(i_new, j_new, k_new);
            assert(new_cid != node.cells[c]->cid);
            g.cells[new_cid].particles.push_back(p);
            p.active = false;
          }
          // erase all inactive particles from current cell
          std::erase_if(node.cells[c]->particles, [](const auto& p) { return !p.active; });
          // todo: could add update to node.active somewhere in here instead of calling g.update_tree() in the push
        }
      } else {
        update_cells(node.children[c], g);
      }
    } // end for(c)
  } // end BorisPush::update_cells()
} // end namespace tf::particles