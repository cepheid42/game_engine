#ifndef GAME_ENGINE_BEAM_SOURCES_HPP
#define GAME_ENGINE_BEAM_SOURCES_HPP

#include "particles.hpp"
#include "vec3.hpp"

#include <adios2.h>

#include <algorithm>

namespace tf::particles
{

struct ParticleBeam {
   explicit ParticleBeam(auto& group_)
   : group(group_)
   {}

   ParticleGroup& group;
   std::vector<Particle> base;

   void apply() {
      group.particles.insert(group.particles.end(), base.begin(), base.end());
      group.sort_particles();
   }
};

auto loadParticleBeam(const auto& filename, auto& group) -> ParticleBeam {
   constexpr vec3 deltas{dx, dy, dz};
   constexpr vec3 mins{x_range[0], y_range[0], z_range[0]};

   std::print("Loading particle source: {}... ", filename);

   ParticleBeam beam(group);

   adios2::ADIOS adios;
   adios2::IO io = adios.DeclareIO("BPReader");
   adios2::Engine reader = io.Open(std::string{sim_path} + filename, adios2::Mode::Read);

   reader.BeginStep();

   const auto p_data = io.InquireVariable<double>("Position");
   const auto v_data = io.InquireVariable<double>("Velocity");
   const auto w_data = io.InquireVariable<double>("Weight");
   const auto g_data = io.InquireVariable<double>("Gamma");

   const auto num_particles = p_data.Shape()[0];

   std::vector<double> p_vec(3 * num_particles);
   std::vector<double> v_vec(3 * num_particles);
   std::vector<double> w_vec(num_particles);
   std::vector<double> g_vec(num_particles);

   reader.Get(p_data, p_vec, adios2::Mode::Sync);
   reader.Get(v_data, v_vec, adios2::Mode::Sync);
   reader.Get(w_data, w_vec, adios2::Mode::Sync);
   reader.Get(g_data, g_vec, adios2::Mode::Sync);

   for (auto i = 0lu; i < num_particles; i++) {
      const vec3 pos{p_vec[3 * i], p_vec[3 * i + 1], p_vec[3 * i + 2]};
      const vec3 vel{v_vec[3 * i], v_vec[3 * i + 1], v_vec[3 * i + 2]};
      const auto weight = w_vec[i];

      const auto loc = ((pos - mins) / deltas);
      const auto gamma = g_vec[i];

      beam.base.emplace_back(
         vel,
         gamma,
         loc,
         loc,
         weight
      );
   }
   reader.EndStep();
   reader.Close();

   std::println("Done.");

   std::ranges::sort(
      beam.base,
      [](const Particle& a, const Particle& b) {
         return morton_encode(getCellIndices<std::size_t>(a.location))
              < morton_encode(getCellIndices<std::size_t>(b.location));
      }
   );

   return beam;
} // end loadParticleBeam()
} // end namespace tf::particles


#endif //GAME_ENGINE_BEAM_SOURCES_HPP