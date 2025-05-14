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
  static std::array<double, 2> shapesCIC(const auto x) {
    return {1.0 - x, x};
  }

  static std::array<double, 3> shapesQuad(const auto x) {
    return {
      0.5 * math::SQR(0.5 - x),
      0.75 - math::SQR(x),
      0.5 * math::SQR(0.5 + x)
    };
  }

  template<std::size_t nx, std::size_t ny, std::size_t nz>
  static std::array<double, nx * ny * nz> computeTotalShape(const auto& xs, const auto& ys, const auto& zs) {
    std::array<double, nx * ny * nz> shapes{};
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
  static auto getField(const auto& emdata, const auto& idxs, const std::array<int, 6>& offsets) {
    std::array<double, nx * ny * nz> fields{};
    int idx = 0;
    for (int i = offsets[0]; i < offsets[1]; i++) {
      for (int j = offsets[2]; j < offsets[3]; j++) {
        // todo: should only run once for 2d, won't work for 3d
        for (int k = offsets[4]; k < offsets[5]; k++) {
          fields[idx] = emdata.template get<F>(idxs[0] + i, idxs[1] + j, idxs[2] + k);
          idx++;
        }
      }
    }
    return fields;
  }

  template<std::size_t N>
  static auto calculateShape(const auto& field, const auto& shape) -> double {
    double result = 0.0;
    for (std::size_t i = 0; i < N; i++) {
      result += field[i] * shape[i];
    }
    return result;
  }

  static std::array<double, 6> FieldAtParticle(Particle& p, const auto& emdata) {
    const auto xr_shapes = shapesCIC(p.location[0]);
    const auto yr_shapes = shapesCIC(p.location[1]);
    const auto zr_shapes = shapesCIC(p.location[2]);
    const auto x_shapes = shapesQuad(p.location[0]);
    // const auto y_shapes = shapesQuad(p.location[1]);
    const auto z_shapes = shapesQuad(p.location[2]);
    constexpr std::array noy = { 1.0 };

    // todo: these shapes are for 2d only, for 3d the 1's need to be replaced
    const auto ex_shapes = computeTotalShape<2, 2, 3>(xr_shapes, yr_shapes, z_shapes);
    const auto ey_shapes = computeTotalShape<3, 1, 3>(x_shapes, noy, z_shapes);
    const auto ez_shapes = computeTotalShape<3, 2, 2>(x_shapes, yr_shapes, zr_shapes);
    const auto bx_shapes = computeTotalShape<2, 1, 3>(xr_shapes, noy, z_shapes);
    const auto by_shapes = computeTotalShape<3, 2, 3>(x_shapes, yr_shapes, z_shapes);
    const auto bz_shapes = computeTotalShape<3, 1, 2>(x_shapes, noy, zr_shapes);

    const auto idxs = morton_decode(p.code);
    const auto Ex_c = getField<0, 2, 2, 3>(emdata, idxs, {0, 2, 0, 2, -1, 2});
    const auto Ey_c = getField<1, 3, 1, 3>(emdata, idxs, {-1, 2, 0, 1, -1, 2});
    const auto Ez_c = getField<2, 3, 2, 2>(emdata, idxs, {-1, 2, 0, 2, 0, 2});

    const auto Bx_c = getField<3, 2, 1, 3>(emdata, idxs, {0, 2, 0, 1, -1, 2});
    const auto By_c = getField<4, 3, 2, 3>(emdata, idxs, {-1, 2, 0, 2, -1, 2});
    const auto Bz_c = getField<5, 3, 1, 2>(emdata, idxs, {-1, 2, 0, 1, 0, 2});

    return {
      calculateShape<12>(Ex_c, ex_shapes),
      calculateShape<9>(Ey_c, ey_shapes),
      calculateShape<12>(Ez_c, ez_shapes),
      calculateShape<6>(Bx_c, bx_shapes),
      calculateShape<18>(By_c, by_shapes),
      calculateShape<6>(Bz_c, bz_shapes)
    };
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
