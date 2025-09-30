#ifndef METRICS_HPP
#define METRICS_HPP

#include "program_params.hpp"
#include "array.hpp"
#include "particles.hpp"
#include "interpolation.hpp"

#include <adios2.h>

#include <unordered_map>
#include <memory>
#include <ctime>

namespace tf::metrics {
// =======================================
// ======== Metrics Base Class ===========
namespace detail {
   struct MetricBase {
      virtual ~MetricBase() = default;
      virtual void write(const std::string&, const std::string&, std::size_t) = 0;
   };
} // end namespace detail

// =====================================
// ======== EM Fields Metric ===========
struct EMFieldsMetric final : detail::MetricBase {
   using pointer_t = Array3D<double>*;
   using field_map = std::unordered_map<std::string, pointer_t>;

   struct FieldVariable {
      pointer_t                   field;
      adios2::Variable<double> variable;
   };

   EMFieldsMetric(const field_map& fields_, adios2::IO&& io_)
   : io(io_)
   {
      for (const auto& [name, field]: fields_) {
         fields.push_back({
            .field = field,
            .variable = io.DefineVariable<double>(
               name,
               {field->nx(), field->ny(), field->nz()}, // shape (global)
               {0, 0, 0},                               // start (local)
               {field->nx(), field->ny(), field->nz()}, // count (local)
               adios2::ConstantDims
            )
         });
      }
   }

   void write(const std::string& dir, const std::string& step_ext, const std::size_t step) override {
      const std::string file{dir + "/fields_" + step_ext};
      static constexpr auto cell_volume = dx * dy * dz;

      io.DefineAttribute<std::size_t>("step", step, "", "/", true);
      io.DefineAttribute<double>("dt", dt);
      io.DefineAttribute<double>("cell_volume", cell_volume);

      adios2::Engine    writer = io.Open(file, adios2::Mode::Write);
      writer.BeginStep();

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
        io.DefineVariable<double>("Position", {group->num_particles(), 3}, {0, 0}, {group->num_particles(), 3})),
     var_vel(
        io.DefineVariable<double>("Velocity", {group->num_particles(), 3}, {0, 0}, {group->num_particles(), 3})),
     var_w(io.DefineVariable<double>("Weight", {group->num_particles(), 1}, {0, 0}, {group->num_particles(), 1})),
     var_gamma(io.DefineVariable<double>("Gamma", {group->num_particles(), 1}, {0, 0}, {group->num_particles(), 1}))
   {
      const auto& nParticles = group->num_particles();
      position.reserve(3 * nParticles);
      velocity.reserve(3 * nParticles);
      weight.reserve(nParticles);
      gamma.reserve(nParticles);
   }

   void write(const std::string& dir, const std::string& step_ext, const std::size_t step) override {
      const std::string file{dir + "/" + group->name + "_dump_" + step_ext};
      const auto& nParticles = group->num_particles();

      io.DefineAttribute<std::string>("name", group->name);
      io.DefineAttribute<std::size_t>("step", step, "", "/", true);
      io.DefineAttribute<double>("dt", dt);
      io.DefineAttribute<double>("mass", group->mass);
      io.DefineAttribute<double>("charge", group->charge);
      // io.DefineAttribute<std::size_t>("atomic_number", group->atomic_number);
      io.DefineAttribute<std::size_t>("num_particles", nParticles, "", "/", true);


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

      static constexpr vec3 delta{dx, dy, dz};
      static constexpr vec3 lb{x_range[0], y_range[0], z_range[0]};

      position.reserve(3 * nParticles);
      velocity.reserve(3 * nParticles);
      weight.reserve(nParticles);
      gamma.reserve(nParticles);

      for (const auto& p: group->particles) {
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

      position.clear();
      velocity.clear();
      weight.clear();
      gamma.clear();
   }

   adios2::IO               io;
   const group_t*           group;
   adios2::Variable<double> var_loc;
   adios2::Variable<double> var_vel;
   adios2::Variable<double> var_w;
   adios2::Variable<double> var_gamma;
   std::vector<double>      position{};
   std::vector<double>      velocity{};
   std::vector<double>      weight{};
   std::vector<double>      gamma{};
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
     var_density(io.DefineVariable<double>("Density", {ncx, ncy, ncz}, {0, 0, 0}, {ncx, ncy, ncz}, adios2::ConstantDims)),
     var_temp(io.DefineVariable<double>("Temperature", {ncx, ncy, ncz}, {0, 0, 0}, {ncx, ncy, ncz}, adios2::ConstantDims)),
     density(ncx, ncy, ncz),
     T_avg(ncx * ncy * ncz),
     KE_total(ncx * ncy * ncz)
   {}

   void update_metrics() {
      using Shape = interp::InterpolationShape<1>::Type;
      static constexpr auto V_cell_inv = 1.0 / (dx * dy * dz);
      static constexpr auto temp_coef  = 2.0 / (3.0 * constants::q_e<double>);

      const auto mc2 = group->mass * constants::c_sqr<double>;

      density.fill(0.0);
      // std::ranges::fill(density, 0.0);
      std::ranges::fill(T_avg, 0.0);
      std::ranges::fill(KE_total, 0.0);


      #pragma omp parallel num_threads(nThreads)
      {
         #pragma omp for
         for (std::size_t pid = 0; pid < group->num_particles(); pid++) {
            const auto& p   = group->particles[pid];
            if (p.disabled) { continue; }
            const vec3 loc_half = particles::getCellIndices<double>(p.location + 1.0) + 0.5;
            const vec3 hid = loc_half.as_type<std::size_t>();

            const vec3 p_half = p.location - loc_half;
            const auto shapeI0 = Shape::shape_array(p_half[0]);
            const auto shapeJ0 = Shape::shape_array(p_half[1]);
            const auto shapeK0 = Shape::shape_array(p_half[2]);

            for (int i = Shape::Begin; i <= Shape::End; ++i) {
               const auto& s0i = shapeI0[i - Shape::Begin];
               for (int j = Shape::Begin; j <= Shape::End; ++j) {
                  const auto& s0j = shapeJ0[j - Shape::Begin];
                  for (int k = Shape::Begin; k <= Shape::End; ++k) {
                     const auto& s0k = shapeK0[k - Shape::Begin];
                     #pragma omp atomic update
                     density(hid[0] + i, hid[1] + j, hid[2] + k) += s0i * s0j * s0k * p.weight;
                  } // end for(k)
               } // end for(j)
            } // end for(i)
         }


         // #pragma omp for
         // for (std::size_t pid = 0; pid < group->num_particles(); pid++) {
         //    const auto& p   = group->particles[pid];
         //    const auto  cid = particles::getCellIndex(p.location);
         //    #pragma omp atomic update
         //    density[cid] += p.weight;
         //    #pragma omp atomic update
         //    KE_total[cid] += p.weight * mc2 * (p.gamma - 1.0);
         // }
         //
         // #pragma omp for
         // for (std::size_t i = 0; i < T_avg.size(); i++) {
         //    if (density[i] == 0.0) { continue; }
         //    T_avg[i] = temp_coef * KE_total[i] / density[i];
         // }
         //
         #pragma omp for simd
         for (std::size_t i = 0; i < density.size(); i++) {
            density[i] *= V_cell_inv;
         }
      } // end parallel
   }

   void write(const std::string& dir, const std::string& step_ext, const std::size_t step) override {
      const std::string file{dir + "/" + group->name + "_" + step_ext};

      update_metrics();

      io.DefineAttribute<std::string>("name", group->name);
      io.DefineAttribute<std::size_t>("step", step, "", "/", true);
      io.DefineAttribute<double>("dt", dt);
      io.DefineAttribute<double>("mass", group->mass);
      io.DefineAttribute<double>("charge", group->charge);
      // io.DefineAttribute<std::size_t>("atomic_number", group->atomic_number);
      io.DefineAttribute<std::size_t>("dims", dims.data(), 3);
      io.DefineAttribute<double>("x_range", x_range.data(), 2);
      io.DefineAttribute<double>("y_range", y_range.data(), 2);
      io.DefineAttribute<double>("z_range", z_range.data(), 2);

      adios2::Engine writer = io.Open(file, adios2::Mode::Write);
      writer.BeginStep();

      writer.Put(var_density, density.data());
      writer.Put(var_temp, T_avg.data());

      writer.EndStep();
      writer.Close();
   } // end write()

   adios2::IO               io;
   const group_t*           group;
   adios2::Variable<double> var_density;
   adios2::Variable<double> var_temp;
   Array3D<double>          density;
   std::vector<double>      T_avg;
   std::vector<double>      KE_total;
   std::array<std::size_t, 3> dims{Nx - 1, Ny - 1, Nz - 1};
};

// =======================================
// ======== Metrics Superclass ===========
class Metrics {
   using metrics_vec = std::vector<std::unique_ptr<detail::MetricBase>>;

public:
   explicit Metrics(std::string data_dir_)
   : data_dir(std::move(data_dir_))
   {
      // const std::time_t time = std::time({});
      // char timestr[std::size("YYMMDDHHMM")];
      // std::strftime(std::data(timestr), std::size(timestr), "%y%m%d%H%M", std::localtime(&time));
      //
      // data_dir += "_" + std::string(timestr);
   }

   void addMetric(std::unique_ptr<detail::MetricBase>&& m) { metrics.push_back(std::move(m)); }

   void write(const std::size_t step) const {
      static constexpr size_t padding      = 10;
      std::string             count_padded = std::to_string(step);
      count_padded.insert(count_padded.begin(), padding - count_padded.length(), '0');
      count_padded += file_ext; // default extension is .bp

      for (const auto& m: metrics) {
         m->write(data_dir, count_padded, step);
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
