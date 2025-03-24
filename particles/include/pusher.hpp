#ifndef PUSHER_HPP
#define PUSHER_HPP

#include "particle.hpp"
#include "em_data.hpp"

#include <cmath>

namespace tf::particles {
  template<typename T>
  T trilinterp(vec3<T>& p, vec3<T>& lo, vec3<T>& hi, std::array<T, 8>& F) {
    const auto tx = (p[0] - lo[0]) / (hi[0] - lo[0]);
    const auto ty = (p[1] - lo[1]) / (hi[1] - lo[1]);
    const auto tz = (p[2] - lo[2]) / (hi[2] - lo[2]);

    // Interpolate in X
    const auto c00 = std::lerp(F[0], F[1], tx);
    const auto c01 = std::lerp(F[2], F[3], tx);
    const auto c10 = std::lerp(F[4], F[5], tx);
    const auto c11 = std::lerp(F[6], F[7], tx);

    // Interpolate in Y
    const auto c0 = std::lerp(c00, c10, ty);
    const auto c1 = std::lerp(c01, c11, ty);

    // Interpolate in Z
    return std::lerp(c0, c1, tz);
  }

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
