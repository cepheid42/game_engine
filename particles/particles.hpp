#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "constants.hpp"
#include "morton.hpp"
#include "program_params.hpp"
#include "vec3.hpp"

#include <boost/sort/sort.hpp>
#include <adios2.h>

#include <span>
#include <utility>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <print>

namespace tf::particles
{
struct Particle {
   vec3<double> location;
   vec3<double> old_location;
   vec3<double> velocity;
   double weight;
   double gamma;
   bool disabled{false};
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
   const auto x = static_cast<std::size_t>(std::floor(loc.x));
   const auto y = static_cast<std::size_t>(std::floor(loc.y));
   const auto z = static_cast<std::size_t>(std::floor(loc.z));
   return z + ((Nz - 1) * y) + ((Ny - 1) * (Nz - 1) * x);
}

constexpr auto calculateGammaV(const auto& v) {
   // Calculates gamma using regular velocity
   return 1.0 / std::sqrt(1.0 - v.length_squared() * constants::over_c_sqr<double>);
}

constexpr auto calculateGammaP(const auto& p) {
   // Calculates gamma using gamma*v (e.g. relativistic momentum but with mass terms canceled)
   return std::sqrt(1.0 + p.length_squared() * constants::over_c_sqr<double>);
}

struct ParticleGroup {
   static constexpr std::size_t SORT_INTERVAL = 100;
   std::string name;
   double mass;
   double charge;
   double qdt_over_2m;
   double initial_y_position{};
   std::vector<Particle> particles{};
   std::map<std::size_t, std::span<Particle>> cell_map{};
   bool cell_map_updated{false};

   ParticleGroup(std::string name_, const double mass_, const double charge_)
      : name(std::move(name_)),
        mass(mass_),
        charge(charge_),
        qdt_over_2m(0.5 * charge * dt / mass) {
   }

   [[nodiscard]] std::size_t num_particles() const { return particles.size(); }

   void update_cell_map() {
      if (cell_map_updated) { return; }
      cell_map.clear();
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


   void reset_y_positions() {
      #pragma omp parallel for simd num_threads(nThreads)
      for (std::size_t pid = 0; pid < particles.size(); pid++) {
         particles[pid].location[1] = initial_y_position;
      }
   }

   void sort_particles() {
      std::erase_if(particles, [](const Particle& p) { return p.disabled; });
      boost::sort::block_indirect_sort(
         particles.begin(), particles.end(),
         [](const Particle& a, const Particle& b)
         {
            return morton_encode(getCellIndices<std::size_t>(a.location)) < morton_encode(
               getCellIndices<std::size_t>(b.location));
         },
         nThreads
      );
      // boost::sort::parallel_stable_sort(
      //    particles.begin(), particles.end(),
      //    [](const Particle& a, const Particle& b) {
      //       return morton_encode(getCellIndices<std::size_t>(a.location)) < morton_encode(getCellIndices<std::size_t>(b.location));
      //    },
      //    nThreads
      // );
   }
}; // end struct ParticleGroup

struct ParticleInitializer {
   static void read_dat(const std::string& filename, auto& group_vec) {
      std::ifstream file(filename);

      if (!file.is_open()) {
         throw std::runtime_error("Particle initialization from file failed: " + filename);
      }

      constexpr vec3 deltas{dx, dy, dz};
      constexpr vec3 mins{x_range[0], y_range[0], z_range[0]};

      std::string comment{};
      std::string name{};
      double mass{};
      double charge{};
      vec3<double> location{};
      vec3<double> velocity{};
      double weight = 0.0;

      std::string line;
      std::getline(file, line);
      std::istringstream buff(line);

      buff >> comment >> name >> mass >> charge;

      ParticleGroup g(name, mass, charge);
      while (std::getline(file, line)) {
         std::istringstream buffer(line);
         buffer >> location >> velocity >> weight;

         location = (location - mins) / deltas;

         // compute Lorentz factor and relativistic momentum
         const auto gamma = calculateGammaV(velocity);

         // add particle to group
         g.particles.emplace_back(
            location,
            location,
            velocity,
            weight,
            gamma,
            false
         );
      }
      file.close();
      if (g.particles.empty()) {
         throw std::runtime_error("Particle initialization failed: Particles vector is empty.");
      }
      g.sort_particles();
      g.initial_y_position = g.particles.empty() ? 0.0 : g.particles[0].location[1];
      group_vec.push_back(g);
   } // end initializeFromFile

   static void read_adios(const std::string& filename, auto& group_vec) {
      constexpr vec3 deltas{dx, dy, dz};
      constexpr vec3 mins{x_range[0], y_range[0], z_range[0]};

      adios2::ADIOS adios;
      adios2::IO io = adios.DeclareIO("BPReader");
      adios2::Engine reader = io.Open(filename, adios2::Mode::Read);

      reader.BeginStep();

      const auto name = io.InquireAttribute<std::string>("name").Data()[0];
      const auto mass = io.InquireAttribute<double>("mass").Data()[0];
      const auto charge = io.InquireAttribute<double>("charge").Data()[0];
      const auto num_particles = io.InquireAttribute<long int>("num_particles").Data()[0];

      ParticleGroup g(name, mass, charge);

      const auto p_data = io.InquireVariable<double>("particles");
      std::vector<double> p_vec(7 * num_particles);

      reader.Get(p_data, p_vec, adios2::Mode::Sync);

      for (auto i = 0lu; i < 7lu * num_particles; i += 7lu) {
         const vec3 pos{p_vec[i], p_vec[i + 1], p_vec[i + 2]};
         const vec3 vel{p_vec[i + 3], p_vec[i + 4], p_vec[i + 5]};
         const auto weight = p_vec[i + 6];

         const auto loc = (pos - mins) / deltas;
         const auto gamma = calculateGammaV(vel);

         g.particles.emplace_back(loc, loc, vel, weight, gamma, false);
      }
      reader.EndStep();
      reader.Close();

      if (g.particles.empty()) {
         throw std::runtime_error("Particle initialization failed: Particles vector is empty.");
      }
      g.sort_particles();
      g.initial_y_position = g.particles.empty() ? 0.0 : g.particles[0].location[1];
      group_vec.push_back(g);
   }

   static void initializeFromFile(const std::string& filename, auto& group_vec) {
      std::print("Loading particle file: {}... ", filename);
      if (filename.ends_with(".bp")) {
         read_adios(filename, group_vec);
      }
      else {
         read_dat(filename, group_vec);
      }
      std::println("Done.");
   }
}; // end struct ParticleInitializer
} // end namespace tf::particles

#endif //PARTICLE_HPP
