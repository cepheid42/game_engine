#ifndef PUSHER_HPP
#define PUSHER_HPP

#include "particle.hpp"
#include "em_data.hpp"

namespace tf::particles {

  template<typename T>
  T linterp(T x, T x0, T x1, T y0, T y1) {
    const auto den = 1.0 / (x1 - x0);
    return den * (y0 * (x1 - x) + y1 * (x - x0));
  }

  template<typename T>
  T bilinterp() {

  }

  compute_t trilinterp() {}

  vec3<compute_t> EFieldAtParticle(Particle& loc, const electromagnetics::EMData& emdata) {
    // todo: trilinear interp via, bilinear interp, via linear interp
  }


  struct BorisIntegrator {
    static void operator()(ParticleGroup& g, electromagnetics::EMData& emdata) {
      // todo: this is where querying an octree would be handy
      for (auto& cell: g.cells) {
        for (auto& chunk : cell) {
          if (chunk.num_active() == 0) { continue; }
          // todo: simd on this loop maybe?
          for (std::size_t pid = 0; pid < chunk.num_active(); pid++) {
            auto& p = chunk.particles[pid];

            const auto eps = g.qdt_over_2m * EFieldAtParticle(p, emdata);
            const auto bet = g.qdt_over_2m * BFieldAtParticle(p, emdata);

            const auto um = p.momentum * g.inv_mass + eps;
            const auto t = bet / std::sqrt(1.0f + um.length_squared() * constants::over_c_sqr);
            const auto s = 2.0f * t / (1.0f + t.length_squared());

            eps += um + cross(um + cross(um, t), s);

            const auto eps_double = static_cast<double>(eps.length_squared());

            p.momentum = g.mass * eps;
            p.gamma = std::sqrt(1.0 + eps_double * constants::over_c_sqr);
          }
        }
      }
    } // end operator()
  }; // end struct BorisIntegrator
} // end namespace tf::particles



#endif //PUSHER_HPP
