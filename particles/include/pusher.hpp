#ifndef PUSHER_HPP
#define PUSHER_HPP

#include "particle.hpp"
#include "em_data.hpp"

#include <cmath>

namespace tf::particles {
  template<typename T>
  T bilinterp(T tx, T ty, const std::array<T, 4>& F) {
    const auto c0 = std::lerp(F[0], F[1], tx);
    const auto c1 = std::lerp(F[2], F[3], ty);

    return std::lerp(c0, c1, ty);
  }

  vec3<compute_t> BFieldAtParticle(const Particle& p,
                                   const vec3<compute_t>& lo,
                                   const vec3<compute_t>& hi,
                                   const std::array<compute_t, 2>& Hx_c,
                                   const std::array<compute_t, 2>& Hy_c,
                                   const std::array<compute_t, 2>& Hz_c)
  {
    const auto t = (p.location - lo) / (hi - lo);
    return {std::lerp(Hx_c[0], Hx_c[1], t[0]),
            std::lerp(Hy_c[0], Hy_c[1], t[1]),
            std::lerp(Hz_c[0], Hz_c[1], t[2])};
  }

  vec3<compute_t> EFieldAtParticle(const Particle& p,
                                   const vec3<compute_t>& lo,
                                   const vec3<compute_t>& hi,
                                   const std::array<compute_t, 4>& Ex_c,
                                   const std::array<compute_t, 4>& Ey_c,
                                   const std::array<compute_t, 4>& Ez_c)
  {
    const auto t = (p.location - lo) / (hi - lo);
    return {bilinterp(t[1], t[2], Ex_c),
            bilinterp(t[0], t[2], Ey_c),
            bilinterp(t[0], t[1], Ez_c)};
  }

  struct BorisPush {
    using emdata_t = electromagnetics::EMData;
    using group_t = ParticleGroup;
    using p_vec = std::vector<Particle>;
    using p_tree = Octree<p_vec>;

    static void update_particle(Particle& p,
                                const ParticleGroup& g,
                                const vec3<compute_t>& lo,
                                const vec3<compute_t>& hi,
                                const std::array<compute_t, 4>& Ex_c,
                                const std::array<compute_t, 4>& Ey_c,
                                const std::array<compute_t, 4>& Ez_c,
                                const std::array<compute_t, 2>& Bx_c,
                                const std::array<compute_t, 2>& By_c,
                                const std::array<compute_t, 2>& Bz_c)
    {
      auto eps = g.qdt_over_2m * EFieldAtParticle(p, lo, hi, Ex_c, Ey_c, Ez_c);
      const auto bet = g.qdt_over_2m * BFieldAtParticle(p, lo, hi, Bx_c, By_c, Bz_c);

      const auto um = p.momentum * g.inv_mass + eps;
      const auto t = bet / static_cast<compute_t>(std::sqrt(1.0 + um.length_squared() * constants::over_c_sqr));
      const auto s = 2.0f * t / (1.0f + t.length_squared());

      eps += um + cross(um + cross(um, t), s);

      const auto eps_double = static_cast<double>(eps.length_squared());

      p.momentum = g.mass * eps;
      p.gamma = std::sqrt(1.0 + eps_double * constants::over_c_sqr);
      p.old_location = p.location;
      p.location += dt * g.inv_mass * p.momentum / static_cast<compute_t>(p.gamma);
    }

    static void update_cell(p_vec& particles, const std::array<std::size_t, 3>& coords, const group_t& g, const emdata_t& emdata) {
      const auto& [i, j, k] = coords;

      const vec3<compute_t> lo{static_cast<compute_t>(i) * dx,
                               static_cast<compute_t>(j) * dy,
                               static_cast<compute_t>(k) * dz};
      const vec3<compute_t> hi{lo[0] + dx, lo[1] + dy, lo[2] + dz};

      const std::array Ex_c = {emdata.Ex(i, j, k), emdata.Ex(i, j + 1, k), emdata.Ex(i, j, k + 1), emdata.Ex(i, j + 1, k + 1)};
      const std::array Ey_c = {emdata.Ey(i, j, k), emdata.Ey(i + 1, j, k), emdata.Ey(i, j, k + 1), emdata.Ey(i + 1, j, k + 1)};
      const std::array Ez_c = {emdata.Ez(i, j, k), emdata.Ez(i, j + 1, k), emdata.Ez(i + 1, j, k), emdata.Ez(i + 1, j + 1, k)};

      const std::array Bx_c = {emdata.Bx(i, j, k), emdata.Bx(i + 1, j, k)};
      const std::array By_c = {emdata.By(i, j, k), emdata.By(i, j + 1, k)};
      const std::array Bz_c = {emdata.Bz(i, j, k), emdata.Bz(i, j, k + 1)};

      for (std::size_t cid = 0; cid < particles.size(); cid++) {
        update_particle(particles[cid], g, lo, hi, Ex_c, Ey_c, Ez_c, Bx_c, By_c, Bz_c);
      }
    }

    static void visit(const p_tree& node, const group_t& g, const emdata_t& emdata) {
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

    static void operator()(const group_t& g, const emdata_t& emdata) {
      // visit(g.tree, g, emdata);
    }
  }; // end struct MomentumIntegrator
} // end namespace tf::particles



#endif //PUSHER_HPP
