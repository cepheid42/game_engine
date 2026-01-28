#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "constants.hpp"
#include "math_utils.hpp"
#include "morton.hpp"
#include "program_params.hpp"
#include "vec3.hpp"

#include <boost/sort/sort.hpp>
#include <adios2.h>

#include <algorithm>
#include <fstream>
#include <print>
#include <span>
#include <utility>
#include <vector>

namespace tf::particles
{
struct Particle {
   vec3<double> velocity; // change to beta and make it a float
   double gamma;
   vec3<double> location;
   vec3<double> old_location;
   double weight;

   [[nodiscard]] bool is_disabled() const { return weight <= 0.0; }
}; // end struct Particle

template <typename T = std::size_t>
constexpr vec3<T> getCellIndices(const auto& loc) {
   return {
      static_cast<T>(std::floor(loc[0])),
      static_cast<T>(std::floor(loc[1])),
      static_cast<T>(std::floor(loc[2]))
   };
}

constexpr std::size_t getCellIndex(const auto& loc) {
   // todo: can probably replace x,y,z with the getCellIndices function above
   // const auto [x, y, z] = getCellIndices(loc);
   const auto x = static_cast<std::size_t>(std::floor(loc.x));
   const auto y = static_cast<std::size_t>(std::floor(loc.y));
   const auto z = static_cast<std::size_t>(std::floor(loc.z));
   return z + ((Nz - 1) * y) + ((Ny - 1) * (Nz - 1) * x);
}

constexpr auto calculateGammaV(const auto& v) {
   // Calculates gamma using regular velocity
   return 1.0 / std::sqrt(1.0 - v.length_squared() * constants::over_c_sqr);
}

constexpr auto calculateGammaP(const auto& p, const auto m) {
   // Calculates gamma using momentum
   return std::sqrt(1.0 + p.length_squared() * constants::over_c_sqr / math::SQR(m));
}

static void initializeFromFile(const std::string& filename, auto& group) {
   constexpr vec3 deltas{dx, dy, dz};
   constexpr vec3 mins{x_range[0], y_range[0], z_range[0]};

   std::print("Loading particle file: {}... ", filename);

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

      group.particles.emplace_back(
         vel,
         gamma,
         loc,
         loc,
         weight
      );
   }
   reader.EndStep();
   reader.Close();
   group.sort_particles();
   std::println("Done.");
}

struct ParticleGroup {
   std::string name;
   std::size_t atomic_number;
   double mass;
   double charge;
   double qdt_over_2m;
   bool cell_map_updated{false};
   bool is_sorted{false};
   std::vector<Particle> particles{};
   std::map<std::size_t, std::span<Particle>> cell_map{};

   ParticleGroup(std::string name_, const std::string& file_, const double mass_, const double charge_, const std::size_t atomic_number_)
   : name(std::move(name_)),
     atomic_number(atomic_number_),
     mass(mass_),
     charge(charge_ * constants::q_e),
     qdt_over_2m((mass != 0.0) ? charge * dt / (2.0 * mass) : 0.0) // for photon groups
   {
      if (not file_.empty()) {
         initializeFromFile(file_, *this);
      }
   }

   explicit ParticleGroup(const ParticleGroupSpec& spec)
   : ParticleGroup(std::string{spec.name}, std::string{spec.filepath}, spec.mass, spec.charge, spec.atomic_number)
   {}

   [[nodiscard]] std::size_t num_particles() const { return particles.size(); }

   void update_cell_map() {
      if (cell_map_updated) { return; }
      cell_map.clear();
      sort_particles();
      auto prev_iter = particles.begin();
      auto prev_code = morton_encode(getCellIndices(prev_iter->location)); // first particles code
      // loop over remaining particles
      for (auto it = particles.begin(); it != particles.end(); ++it) {
         const auto cur_code = morton_encode(getCellIndices(it->location));
         if (prev_code != cur_code) {
            cell_map.insert({prev_code, std::span{prev_iter, it}});
            prev_code = cur_code;
            prev_iter = it;
         }
      }
      cell_map.insert_or_assign(prev_code, std::span{prev_iter, particles.end()});
      cell_map_updated = true;
   }

   void reset_positions() {
      if constexpr (x_collapsed or y_collapsed or z_collapsed) {
         #pragma omp parallel for simd num_threads(nThreads)
         for (std::size_t pid = 0; pid < particles.size(); pid++) {
            if constexpr (x_collapsed) { particles[pid].location[0] = 0.5; }
            if constexpr (y_collapsed) { particles[pid].location[1] = 0.5; }
            if constexpr (z_collapsed) { particles[pid].location[2] = 0.5; }
         }
      }
   }

   void sort_particles() {
      if (is_sorted) { return; }
      std::erase_if(particles, [](const Particle& p) { return p.is_disabled(); });
      boost::sort::block_indirect_sort(
         particles.begin(), particles.end(),
         [](const Particle& a, const Particle& b) {
            return morton_encode(getCellIndices<std::size_t>(a.location))
                 < morton_encode(getCellIndices<std::size_t>(b.location));
         },
         nThreads
      );
      is_sorted = true;
   }
}; // end struct ParticleGroup
} // end namespace tf::particles

#endif //PARTICLE_HPP
