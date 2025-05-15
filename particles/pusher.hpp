#ifndef PUSHER_HPP
#define PUSHER_HPP

#include "program_params.hpp"
#include "em_params.hpp"
#include "particles.hpp"
#include "em_data.hpp"
#include "constants.hpp"
#include "morton.hpp"

#include <cmath>

namespace tf::particles {
  template<int D>
  auto rotateOrigin(const auto& p) {
    if constexpr (D == 0) {
      return decltype(p){p[1], p[2], p[0]};
    }
    else if constexpr (D == 1) {
      return decltype(p){p[2], p[0], p[1]};
    }
    else {
      return p;
    }
  }

  auto shapesQuad(const auto x) {
      const auto absx = std::abs(x);
    if (absx < 0.5) {
      return 0.75 - math::SQR(x);
    }
    return 0.5 * math::SQR(1.5 - absx);
  }

  template<int D>
  auto EFieldToParticle(const auto& E, const auto& p0, const auto& cids) {
    static constexpr int i0 = D == 0 ? 0 : -1;
    static constexpr int i1 = 1;

    static constexpr int j0 = -1;
    static constexpr int j1 = 1;

    static constexpr int k0 = -1;
    static constexpr int k1 = 1;

    // rotate elements to match loop structure
    const auto [ci, cj, ck] = rotateOrigin<D>(cids);
    const auto [x0, y0, z0] = rotateOrigin<D>(p0);

    auto result = 0.0_fp;
    for (int i = i0; i <= i1; ++i) {
      const auto s0i = shape(x0 - i);

      for (int j = j0; j <= j1; ++j) {
        const auto s0j = shape(y0 - j);

        for (int k = k0; k <= k1; ++k) {
          const auto s0k = shape(z0 - k);

          // undo the rotation to get proper indices back
          const auto [x, y, z] = rotateOrigin<D == 2 ? D : !D>(vec3{i + ci, j + cj, k + ck});

          result += s0i * s0j * s0k * E.template get<D>(ci + i, cj + j, ck + k);
        }
      }
    }
    return result;
  }

  template<int D>
  auto BFieldToParticle(const auto& B, const auto& p0, const auto& cids) {
    static constexpr int i0 = D == 0 ? 0 : -1;
    static constexpr int i1 = 1;

    static constexpr int j0 = -1;
    static constexpr int j1 = 1;

    static constexpr int k0 = -1;
    static constexpr int k1 = 1;

    // rotate elements to match loop structure
    const auto [ci, cj, ck] = rotateOrigin<D>(cids);
    const auto [x0, y0, z0] = rotateOrigin<D>(p0);

    auto result = 0.0_fp;
    for (int i = i0; i <= i1; ++i) {
      const auto s0i = shape(x0 - i);

      for (int j = j0; j <= j1; ++j) {
        const auto s0j = shape(y0 - j);

        for (int k = k0; k <= k1; ++k) {
          const auto s0k = shape(z0 - k);

          // undo the rotation to get proper indices back
          const auto [x, y, z] = rotateOrigin<D == 2 ? D : !D>(vec3{i + ci, j + cj, k + ck});

          result += s0i * s0j * s0k * E(ci + i, cj + j, ck + k);
        }
      }
    }
    return result;
  }

  static std::array<double, 6> FieldAtParticle(Particle& p, const auto& emdata) {
    const auto cids = morton_decode(p.code);
    const auto exc = EFieldToParticle<0>(emdata.Ex, p.location, cids);

    const auto bxc = BFieldToParticle<0>(emdata.Bx, p.location, cids);

    return {};
  } // end FieldAtParticle

  struct BorisPush {
    using emdata_t = electromagnetics::EMData;
    using group_t = ParticleGroup;

    static constexpr std::size_t BC_DEPTH = PMLDepth - 5;

    static void update_velocity(Particle& p, const emdata_t& emdata, auto qdt) {
      const auto emf = FieldAtParticle(p, emdata);

      const vec3<double> eps{qdt * emf[0], qdt * emf[1], qdt * emf[2]};
      const vec3<double> bet{qdt * emf[3], qdt * emf[4], qdt * emf[5]};

      const auto um = p.velocity.as_type<double>() + eps;
      const auto t = bet / std::sqrt(1.0 + um.length_squared() * constants::over_c_sqr<double>);
      const auto s = 2.0 * t / (1.0 + t.length_squared());

      const auto v = eps + um + cross(um + cross(um, t), s);
      const auto gamma = std::sqrt(1.0 + v.length_squared() * constants::over_c_sqr<double>);

      p.gamma = gamma;
      p.velocity = v.as_type<compute_t>();
    } // end update_velocity()

    // #pragma omp declare simd notinbranch
    static void update_position(Particle& p) {
      static constexpr vec3 delta_inv{dt / dx, dt / dy, dt / dz};

      const auto [iold, jold, kold] = morton_decode(p.code);
      const auto new_loc = p.location + delta_inv * p.velocity / p.gamma;
      const vec3 offsets = {std::floor(new_loc[0]),
                            std::floor(new_loc[1]),
                            std::floor(new_loc[2])};

      assert(offsets[1] == 0.0);

      const auto inew = iold + static_cast<std::size_t>(offsets[0]);
      const auto jnew = jold + static_cast<std::size_t>(offsets[1]);
      const auto knew = kold + static_cast<std::size_t>(offsets[2]);

      if (inew < BC_DEPTH or inew > Nx - BC_DEPTH or knew < BC_DEPTH or knew > Nz - BC_DEPTH) {
        p.code = group_t::DISABLED;
      } else {
        p.code = morton_encode(inew, jnew, knew);
        std::println("{}, {}, {}", inew, jnew, knew);
      }

      p.old_location = p.location - offsets;
      p.location = new_loc - offsets;
    } // end update_position()


    static void advance_velocity(group_t& g, const emdata_t& emdata) {
      // #pragma omp parallel for num_threads(nThreads)
      for (std::size_t pid = 0; pid < g.num_particles(); pid++) {
        update_velocity(g.particles[pid], emdata, g.qdt_over_2m);
      }
    } // end advance_velocity

    static void advance_position(group_t& g) {
      // #pragma omp parallel for simd num_threads(nThreads)
      for (std::size_t pid = 0; pid < g.num_particles(); pid++) {
        update_position(g.particles[pid]);
      }
      std::erase_if(g.particles, [](const Particle& p) { return p.code == group_t::DISABLED; });
    } // end advance_position


    static void backstep_velocity(group_t& g, const emdata_t& emdata) {
      // Easiest way is to just copy the velocity update and make qdt negative
      // #pragma omp parallel for num_threads(nThreads)
      for (std::size_t pid = 0; pid < g.num_particles(); pid++) {
        update_velocity(g.particles[pid], emdata, -0.5 * g.qdt_over_2m);
      }
    }

    static void operator()(group_t& g, const emdata_t& emdata, const std::size_t step) {
      advance_velocity(g, emdata);
      advance_position(g);

      if (step % group_t::SORT_INTERVAL == 0) {
        g.sort_particles();
      }
    } // end operator()
  }; // end struct BorisPush
} // end namespace tf::particles


#endif //PUSHER_HPP
