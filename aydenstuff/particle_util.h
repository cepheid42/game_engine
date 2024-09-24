//
// Created by akis on 7/29/24.
//

#ifndef PARTICLE_UTIL_H
#define PARTICLE_UTIL_H

#include <concepts>
#include <vector>

#include "vector2d.h"
#include "vector3d.h"

namespace tf_util::particle
{
  template <std::floating_point fp>
  auto computeLorentzFactorP(const vec3<fp>& momentum, const fp mass) {
    return std::sqrt( 1.0 + momentum.length_squared() / math::SQR(mass * constants::c) );
  } // end computeLorentzFactor

  template <std::floating_point fp>
  auto computeLorentzFactorV(const vec3<fp>& velocity) {
    return 1.0 / std::sqrt( 1.0 - velocity.length_squared() / constants::c_sqr);
  } // end computeLorentzFactor

  template<std::integral I, std::floating_point fp>
  struct TrajectorySegment {
    TrajectorySegment(fp weight, vec2<fp> p0, vec2<fp> p1, vec2<I> cid) : weight(weight), p0(p0), p1(p1), cid(cid) {}

    fp weight{};
    vec2<fp> p0{}, p1{};
    vec2<I> cid{};
  };

  template<std::integral I, std::floating_point fp>
  using Trajectory = std::vector<TrajectorySegment<I, fp>>;

  static Trajectory<int> split_trajectory(particle_t& particle, SimulationMesh& mesh)
    {
      auto& reference_mesh = interpolation::referenceMesh<shape_t>(mesh.full_mesh);

      auto stencil0 = interpolation::findStencilReferenceIndex<shape_t>(particle.old_location, mesh.full_mesh);
      auto stencil1 = interpolation::findStencilReferenceIndex<shape_t>(particle.location, mesh.full_mesh);

      vec2<int> offset{
        static_cast<int>(stencil1[0]) - static_cast<int>(stencil0[0]),
        static_cast<int>(stencil1[1]) - static_cast<int>(stencil0[1])
      };

      const auto num_crossings = static_cast<size_t>(offset[0] != 0) + static_cast<size_t>(offset[1] != 0);

      Trajectory<int> trajectory(0);
      if (num_crossings == 0) {
        trajectory.emplace_back(particle.weight, particle.old_location, particle.location, stencil0);
      } else {
        auto t_stack = std::vector<fptype>{1.0};

        if (offset[0] != 0) {
          auto i_crossing = offset[0] < 0 ? stencil0[0] : stencil1[0];
          auto x_crossing = reference_mesh.x[i_crossing];
          t_stack.push_back((x_crossing - particle.old_location[0]) / (particle.location[0] - particle.old_location[0]));
        }

        if (offset[1] != 0) {
          auto k_crossing = offset[1] < 0 ? stencil0[1] : stencil1[1];
          auto z_crossing = reference_mesh.z[k_crossing];
          t_stack.push_back((z_crossing - particle.old_location[1]) / (particle.location[1] - particle.old_location[1]));
        }

        std::ranges::sort(t_stack, std::ranges::greater());
        fptype t_start = 0.0;
        while (!t_stack.empty()) {
          auto t_end = t_stack.back();
          t_stack.pop_back();

          auto p0 = (1.0 - t_start) * particle.old_location + t_start * particle.location;
          auto p1 = (1.0 - t_end) * particle.old_location + t_end * particle.location;

          trajectory.emplace_back(particle.weight, p0, p1, interpolation::findStencilReferenceIndex<shape_t>(p0, mesh.full_mesh));

          t_start = t_end;
        } // end while(t_stack)
        //
      } // endif (num_crosssings)

      return trajectory;
      //
    } // end split_trajectory

  //
} // end tf_util::particle

#endif //PARTICLE_UTIL_H
