#ifndef METRICS_HPP
#define METRICS_HPP

#include "program_params.hpp"
#include "array.hpp"
#include "particles.hpp"

#include <unordered_map>
#include <memory>

#include <adios2.h>

namespace tf::metrics {
// =======================================
// ======== Metrics Base Class ===========
namespace detail {
   struct MetricBase {
      virtual ~MetricBase() = default;
      virtual void write(const std::string&, const std::string&) = 0;
   };
} // end namespace detail

// =====================================
// ======== EM Fields Metric ===========
struct EMFieldsMetric final : detail::MetricBase {
   using pointer_t = Array3D<compute_t>*;
   using field_map = std::unordered_map<std::string, pointer_t>;

   struct FieldVariable {
      pointer_t                   field;
      adios2::Variable<compute_t> variable;
   };

   EMFieldsMetric(const field_map& fields_, adios2::IO&& io_)
   : io(io_)
   {
      for (const auto& [name, field]: fields_) {
         fields.push_back({
            .field = field,
            .variable = io.DefineVariable<compute_t>(
               name,
               {field->nx(), field->ny(), field->nz()}, // shape (global)
               {0, 0, 0},                               // start (local)
               {field->nx(), field->ny(), field->nz()}, // count (local)
               adios2::ConstantDims
            )
         });
      }
   }

   void write(const std::string& dir, const std::string& step_ext) override {
      const std::string file{dir + "/fields_" + step_ext};
      adios2::Engine    writer = io.Open(file, adios2::Mode::Write);
      writer.BeginStep();

      // ReSharper disable once CppUseElementsView
      for (auto& [field, variable]: fields) {
         writer.Put(variable, field->data());
      }

      writer.EndStep();
      writer.Close();
   }

   adios2::IO                 io;
   std::vector<FieldVariable> fields;
};

// =========================================
// ======== Particle Dump Metric ===========
struct ParticleDumpMetric final : detail::MetricBase {
   using group_t = particles::ParticleGroup;

   ParticleDumpMetric(const group_t* group_, adios2::IO&& io_)
   : io(io_),
     group(group_),
     var_loc(
        io.DefineVariable<compute_t>("Position", {group->num_particles(), 3}, {0, 0}, {group->num_particles(), 3})),
     var_vel(
        io.DefineVariable<compute_t>("Velocity", {group->num_particles(), 3}, {0, 0}, {group->num_particles(), 3})),
     var_w(io.DefineVariable<compute_t>("Weight", {group->num_particles(), 1}, {0, 0}, {group->num_particles(), 1})),
     var_gamma(io.DefineVariable<double>("Gamma", {group->num_particles(), 1}, {0, 0}, {group->num_particles(), 1}))
   {}

   void write(const std::string& dir, const std::string& step_ext) override {
      const std::string file{dir + "/" + group->name + "_dump_" + step_ext};

      // io.DefineAttribute<std::string>("name", group->name);
      // io.DefineAttribute<std::size_t>("num_particles", group->num_particles);
      // io.DefineAttribute<std::size_t>("atomic_number", group->atomic_number);
      // io.DefineAttribute<compute_t>("mass", group->mass);
      // io.DefineAttribute<compute_t>("charge", group->charge);
      const auto& nParticles = group->num_particles();

      var_loc.SetShape({nParticles, 3});
      var_loc.SetSelection({{0, 0}, {nParticles, 3}}); // {{start}, {count}}

      var_vel.SetShape({nParticles, 3});
      var_vel.SetSelection({{0, 0}, {nParticles, 3}}); // {{start}, {count}}

      var_w.SetShape({nParticles, 1});
      var_w.SetSelection({{0, 0}, {nParticles, 1}}); // {{start}, {count}}

      var_gamma.SetShape({nParticles, 1});
      var_gamma.SetSelection({{0, 0}, {nParticles, 1}}); // {{start}, {count}}

      adios2::Engine writer = io.Open(file, adios2::Mode::Write);
      writer.BeginStep();

      constexpr vec3 delta{dx, dy, dz};
      constexpr vec3 lb{x_range[0], y_range[0], z_range[0]};

      std::vector<compute_t> position{};
      std::vector<compute_t> velocity{};
      std::vector<compute_t> weight{};
      std::vector<double>    gamma{};

      position.reserve(3 * nParticles);
      velocity.reserve(3 * nParticles);
      weight.reserve(nParticles);
      gamma.reserve(nParticles);

      for (const auto& p: group->particles) {
         // const auto idxs = particles::getCIDs(p.location);
         for (std::size_t d = 0; d < 3; d++) {
            position.push_back(lb[d] + delta[d] * p.location[d]);
            velocity.push_back(p.velocity[d]);
         }
         weight.push_back(p.weight);
         gamma.push_back(p.gamma);
      }

      writer.Put(var_loc, position.data());
      writer.Put(var_vel, velocity.data());
      writer.Put(var_w, weight.data());
      writer.Put(var_gamma, gamma.data());

      writer.EndStep();
      writer.Close();
   }

   adios2::IO                  io;
   const group_t*              group;
   adios2::Variable<compute_t> var_loc;
   adios2::Variable<compute_t> var_vel;
   adios2::Variable<compute_t> var_w;
   adios2::Variable<double>    var_gamma;
};

// ========================================================
// ======== Particle Density/Temperature Metric ===========
struct ParticleMetric final : detail::MetricBase {
   using group_t = particles::ParticleGroup;

   ParticleMetric(const group_t*    g_,
                  adios2::IO&&      io_,
                  const std::size_t ncx,
                  const std::size_t ncy,
                  const std::size_t ncz)
   : io(io_),
     group(g_),
     var_density(io.DefineVariable<compute_t>("Density", {ncx, ncy, ncz}, {0, 0, 0}, {ncx, ncy, ncz},
                                              adios2::ConstantDims)),
     var_temp(io.DefineVariable<compute_t>("Temperature", {ncx, ncy, ncz}, {0, 0, 0}, {ncx, ncy, ncz},
                                           adios2::ConstantDims)),
     density(Ncx * Ncy * Ncz),
     T_avg(Ncx * Ncy * Ncz),
     KE_total(Ncx * Ncy * Ncz)
   {}

   void update_metrics() {
      static constexpr auto V_cell_inv = 1.0_fp / (dx * dy * dz);
      static constexpr auto temp_coef  = 2.0_fp / (3.0_fp * constants::q_e<compute_t>);

      const auto mc2 = group->mass * constants::c_sqr<compute_t>;

      std::ranges::fill(density, 0.0_fp);
      std::ranges::fill(T_avg, 0.0_fp);
      std::ranges::fill(KE_total, 0.0_fp);

      #pragma omp parallel num_threads(nThreads)
      {
         #pragma omp for
         for (std::size_t pid = 0; pid < group->num_particles(); pid++) {
            const auto& p         = group->particles[pid];
            const auto  [i, j, k] = particles::getCIDs(p.location);
            const auto  cid       = get_cid(i, j, k);
            #pragma omp atomic update
            density[cid] += p.weight;
            #pragma omp atomic update
            KE_total[cid] += p.weight * mc2 * (p.gamma - 1.0);
         }

         #pragma omp for
         for (std::size_t i = 0; i < T_avg.size(); i++) {
            if (density[i] == 0.0_fp) { continue; }
            T_avg[i] = temp_coef * KE_total[i] / density[i];
         }

         #pragma omp for
         for (std::size_t i = 0; i < density.size(); i++) {
            density[i] *= V_cell_inv;
         }
      }
   }

   void write(const std::string& dir, const std::string& step_ext) override {
      update_metrics();

      const std::string file{dir + "/" + group->name + "_" + step_ext};
      adios2::Engine    writer = io.Open(file, adios2::Mode::Write);
      writer.BeginStep();

      writer.Put(var_density, density.data());
      writer.Put(var_temp, T_avg.data());

      writer.EndStep();
      writer.Close();
   } // end write()

   adios2::IO                  io;
   const group_t*              group;
   adios2::Variable<compute_t> var_density;
   adios2::Variable<compute_t> var_temp;
   std::vector<compute_t>      density;
   std::vector<compute_t>      T_avg;
   std::vector<compute_t>      KE_total;
};

// =======================================
// ======== Metrics Superclass ===========
class Metrics {
   using metrics_vec = std::vector<std::unique_ptr<detail::MetricBase>>;

public:
   explicit Metrics(std::string data_dir_)
   : data_dir(std::move(data_dir_))
   {}

   void addMetric(std::unique_ptr<detail::MetricBase>&& m) { metrics.push_back(std::move(m)); }

   void write(const std::size_t step) const {
      static constexpr size_t padding      = 10;
      std::string             count_padded = std::to_string(step);
      count_padded.insert(count_padded.begin(), padding - count_padded.length(), '0');
      count_padded += file_ext; // default extension is .bp

      for (const auto& m: metrics) {
         m->write(data_dir, count_padded);
      }
   }

public:
   adios2::ADIOS adios{};

private:
   std::string file_ext{".bp"};
   std::string data_dir{};
   metrics_vec metrics{};
};
}

#endif //METRICS_HPP
