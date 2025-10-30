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
   vec3<double> velocity;
   double gamma;
   vec3<float> location;
   vec3<float> old_location;
   float weight; // does this NEED to be a double?

   [[nodiscard]] bool is_disabled() const { return weight < 0.0f; }
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

// constexpr auto calculateGamma(const double v2) {
//    return 1.0 / std::sqrt(1.0 - v2 * constants::over_c_sqr<double>);
// }

constexpr auto calculateGammaV(const auto& v) {
   // Calculates gamma using regular velocity
   return 1.0 / std::sqrt(1.0 - v.length_squared() * constants::over_c_sqr<double>);
}

constexpr auto calculateGammaP(const auto& p, const auto m) {
   // Calculates gamma using gamma*v (e.g. relativistic momentum but with mass terms canceled)
   return std::sqrt(1.0 + p.length_squared() / math::SQR(m) * constants::over_c_sqr<double>);
}

struct ParticleGroup {
   static constexpr std::size_t SORT_INTERVAL = 50;
   std::string name;
   double mass;
   double charge;
   std::size_t atomic_number;
   double qdt_over_2m;
   vec3<float> initial_position{};
   bool tracer;
   bool sourcer;
   std::vector<Particle> particles{};
   std::map<std::size_t, std::span<Particle>> cell_map{};
   bool cell_map_updated{false};

   ParticleGroup(std::string name_, const double mass_, const double charge_, const std::size_t atomic_number_, const bool tracer_flag = false, const bool sourcer_flag = false)
   : name(std::move(name_)),
     mass(mass_),
     charge(charge_),
     atomic_number(atomic_number_),
     qdt_over_2m(0.5 * charge * dt / mass),
     tracer(tracer_flag),
     sourcer(sourcer_flag)
   {}

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

   static void reset_positions() {}

   void reset_positions() requires (x_collapsed or y_collapsed or z_collapsed) {
      #pragma omp parallel for simd num_threads(nThreads)
      for (std::size_t pid = 0; pid < particles.size(); pid++) {
         if constexpr (x_collapsed) {particles[pid].location[0] = initial_position.x;}
         if constexpr (y_collapsed) {particles[pid].location[1] = initial_position.y;}
         if constexpr (z_collapsed) {particles[pid].location[2] = initial_position.z;}
      }
   }

   void sort_particles() {
      std::erase_if(particles, [](const Particle& p) { return p.is_disabled(); });
      boost::sort::block_indirect_sort(
         particles.begin(), particles.end(),
         [](const Particle& a, const Particle& b) {
            return morton_encode(getCellIndices<std::size_t>(a.location))
                 < morton_encode(getCellIndices<std::size_t>(b.location));
         },
         nThreads
      );
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
      vec3<double> pos{};
      vec3<double> velocity{};
      float weight = 0.0;

      std::string line;
      std::getline(file, line);
      std::istringstream buff(line);

      buff >> comment >> name >> mass >> charge;

      ParticleGroup g(name, mass, charge, 0);
      while (std::getline(file, line)) {
         std::istringstream buffer(line);
         buffer >> pos >> velocity >> weight;

         const auto loc = ((pos - mins) / deltas).as_type<float>();

         // compute Lorentz factor and relativistic momentum
         const auto gamma = calculateGammaV(velocity);

         // add particle to group
         g.particles.emplace_back(
            velocity,
            gamma,
            pos.as_type<float>(),
            pos.as_type<float>(),
            weight
         );
      }
      file.close();
      if (g.particles.empty()) {
         throw std::runtime_error("Particle initialization failed: Particles vector is empty.");
      }
      g.sort_particles();
      g.initial_position = g.particles.empty() ? vec3<float>{} : g.particles[0].location;
      group_vec.push_back(g);
   } // end initializeFromFile

   static void read_adios(const std::string& filename, auto& group_vec) {
      constexpr vec3 deltas{dx, dy, dz};
      constexpr vec3 mins{x_range[0], y_range[0], z_range[0]};

      adios2::ADIOS adios;
      adios2::IO io = adios.DeclareIO("BPReader");
      adios2::Engine reader = io.Open(filename, adios2::Mode::Read);

      reader.BeginStep();

      const auto name = io.InquireAttribute<std::string>("Name").Data()[0];
      const auto mass = io.InquireAttribute<double>("Mass").Data()[0];
      const auto charge = io.InquireAttribute<double>("Charge").Data()[0];
      const std::size_t atomic_number = io.InquireAttribute<unsigned long int>("Atomic Number").Data()[0];
      const bool tracer = static_cast<bool>(io.InquireAttribute<unsigned long int>("Tracer").Data()[0]);
      const bool sourcer = static_cast<bool>(io.InquireAttribute<unsigned long int>("Sourcer").Data()[0]);

      ParticleGroup g(name, mass, charge, atomic_number, tracer, sourcer);

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
         const auto weight = static_cast<float>(w_vec[i]);

         const auto loc = ((pos - mins) / deltas);
         // const auto gamma = calculateGamma(vel);
         const auto gamma = g_vec[i];

         g.particles.emplace_back(
            vel,
            gamma,
            loc.as_type<float>(),
            loc.as_type<float>(),
            weight
         );
      }
      reader.EndStep();
      reader.Close();

      if (g.particles.empty()) {
         throw std::runtime_error("Particle initialization failed: Particles vector is empty.");
      }
      g.sort_particles();
      g.initial_position = g.particles.empty() ? vec3<float>{} : g.particles[0].location;
      group_vec.push_back(g);
   }

   static void initializeFromFile(const std::string& filename, auto& group_vec) {
      std::print("Loading particle file: {}... ", filename);
      if (filename.ends_with(".bp")) {
         read_adios(filename, group_vec);
      } else {
         read_dat(filename, group_vec);
      }
      std::println("Done.");
   }
}; // end struct ParticleInitializer
} // end namespace tf::particles

#endif //PARTICLE_HPP
