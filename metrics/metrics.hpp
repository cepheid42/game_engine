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
      virtual void write(const std::string&, const std::string&, std::size_t, double) = 0;
      std::size_t interval = 100lu;
   };
} // end namespace detail

// =====================================
// ====== EM Field Slice Metric ========
struct EMFieldsMetric final : detail::MetricBase {
   using pointer_t = Array3D<double>*;
   using emdata_t = electromagnetics::EMData;

   struct FieldVariable {
      pointer_t                field;
      adios2::Variable<double> variable;
   };

   EMFieldsMetric(auto& em_map, adios2::IO&& io_)
   : io(io_),
     var_step(io.DefineVariable<std::size_t>("Step")),
     var_dt(io.DefineVariable<double>("dt")),
     var_time(io.DefineVariable<double>("Time"))
     // BFields({
     //    {"Bx", Array3D<double>{em_map["Hx"].dims()}},
     //    {"By", Array3D<double>{em_map["Hy"].dims()}},
     //    {"Bz", Array3D<double>{em_map["Hz"].dims()}}
     // })
   {
      BFields.insert({"Bx", Array3D<double>(Nx, Ny - 1, Nz - 1)});
      BFields.insert({"By", Array3D<double>(Nx - 1, Ny, Nz - 1)});
      BFields.insert({"Bz", Array3D<double>(Nx - 1, Ny - 1, Nz)});

      for (auto& [name, field] : em_map) {
         auto nm = name;
         if (nm[0] == 'H') {
            nm[0] = 'B';
         }
         const auto [xs, ys, zs] = field.dims();
         fields.push_back(FieldVariable{
            .field = &field,
            .variable = io.DefineVariable<double>(
               nm,
               {xs, ys, zs}, // shape
               {0, 0, 0},    // start
               {xs, ys, zs}, // count
               adios2::ConstantDims
            )
         });
      }

      io.DefineAttribute<std::string>("Unit", "s", "Time");
      io.DefineAttribute<std::string>("Unit", "s", "dt");
      io.DefineAttribute<std::string>("Cell Volume/Unit", "m^{3}");
      io.DefineAttribute<double>("Cell Volume", dx * dy * dz);
   }

   static void write(const auto&, const auto&, const auto, const auto) requires(!em_enabled) {}

   void write(const std::string& dir, const std::string& step_ext, const std::size_t step, const double time) override {
      const std::string file = dir + "/fields_" + step_ext;

      if (step % em_save_interval != 0) { return; }

      adios2::Engine writer = io.Open(file, adios2::Mode::Write);
      writer.BeginStep();

      for (auto& [field, var] : fields) {
         if (var.Name()[0] == 'B') {
            std::ranges::transform(field->begin(), field->end(), BFields[var.Name()].begin(), [](const double el){ return constants::mu0<double> * el; });
            field = &BFields[var.Name()];
         }

         writer.Put(var, field->data());
         writer.Put(var_step, step);
         writer.Put(var_dt, dt);
         writer.Put(var_time, time);
      }

      writer.EndStep();
      writer.Close();
   }

   adios2::IO                    io;
   adios2::Variable<std::size_t> var_step;
   adios2::Variable<double>      var_dt;
   adios2::Variable<double>      var_time;
   std::vector<FieldVariable>    fields;
   std::unordered_map<std::string, Array3D<double>> BFields;
}; // end struct EMFieldsMetric

// =====================================
// ======== EM Energy Metric ===========
struct EMTotalEnergyMetric final : detail::MetricBase {
   using data_t = Array3D<double>;
   using field_map = std::unordered_map<std::string, data_t*>;

   EMTotalEnergyMetric(field_map& em_map_, adios2::IO&& io_)
   : io(io_),
     em_map(em_map_),
     var_ex_energy(io.DefineVariable<double>("Ex Energy")),
     var_ey_energy(io.DefineVariable<double>("Ey Energy")),
     var_ez_energy(io.DefineVariable<double>("Ez Energy")),
     var_bx_energy(io.DefineVariable<double>("Bx Energy")),
     var_by_energy(io.DefineVariable<double>("By Energy")),
     var_bz_energy(io.DefineVariable<double>("Bz Energy")),
     var_step(io.DefineVariable<std::size_t>("Step")),
     var_dt(io.DefineVariable<double>("dt")),
     var_time(io.DefineVariable<double>("Time"))
   {
      io.DefineAttribute<std::string>("Unit", "s", "Time");
      io.DefineAttribute<std::string>("Unit", "s", "dt");
      io.DefineAttribute<double>("Cell Volume", dx * dy * dz);
   }

   static void write(const auto&, const auto&, const auto, const auto) requires(!em_enabled) {}

   void write(const std::string& dir, const std::string&, const std::size_t step, const double time) override {
      static constexpr auto cell_volume = dx * dy * dz;
      const std::string file{dir + "/fields_energy.bp"};

      if (step % em_save_interval != 0) { return; }

      auto Ex_energy = 0.0;
      auto Ey_energy = 0.0;
      auto Ez_energy = 0.0;
      auto Bx_energy = 0.0;
      auto By_energy = 0.0;
      auto Bz_energy = 0.0;

      // todo: fix this to use the right volumes for non-uniform mesh!
      #pragma omp parallel num_threads(nThreads) default(shared)
      {
         #pragma omp for reduction(+:Ex_energy)
         for (auto i = 0zu; i < em_map["Ex"]->size(); i++) {
           Ex_energy += math::SQR(em_map["Ex"]->operator[](i)) * cell_volume;
        }
         #pragma omp for reduction(+:Ey_energy)
         for (auto i = 0zu; i < em_map["Ey"]->size(); i++) {
            Ey_energy += math::SQR(em_map["Ey"]->operator[](i)) * cell_volume;
         }
         #pragma omp for reduction(+:Ez_energy)
         for (auto i = 0zu; i < em_map["Ez"]->size(); i++) {
            Ez_energy += math::SQR(em_map["Ez"]->operator[](i)) * cell_volume;
         }
         #pragma omp for reduction(+:Bx_energy)
         for (auto i = 0zu; i < em_map["Hx"]->size(); i++) {
            Bx_energy += math::SQR(em_map["Hx"]->operator[](i)) * cell_volume;
         }
         #pragma omp for reduction(+:By_energy)
         for (auto i = 0zu; i < em_map["Hy"]->size(); i++) {
            By_energy += math::SQR(em_map["Hy"]->operator[](i)) * cell_volume;
         }
         #pragma omp for reduction(+:Bz_energy)
         for (auto i = 0zu; i < em_map["Hz"]->size(); i++) {
            Bz_energy += math::SQR(em_map["Hz"]->operator[](i)) * cell_volume;
         }
      }

      // append after step 0
      adios2::Engine writer = io.Open(file, (step == 0) ? adios2::Mode::Write : adios2::Mode::Append);
      writer.BeginStep();
      writer.Put(var_ex_energy, Ex_energy * 0.5 * tf::constants::eps0<double>);
      writer.Put(var_ey_energy, Ey_energy * 0.5 * tf::constants::eps0<double>);
      writer.Put(var_ez_energy, Ez_energy * 0.5 * tf::constants::eps0<double>);
      writer.Put(var_bx_energy, Bx_energy * 0.5 * tf::constants::mu0<double>);
      writer.Put(var_by_energy, By_energy * 0.5 * tf::constants::mu0<double>);
      writer.Put(var_bz_energy, Bz_energy * 0.5 * tf::constants::mu0<double>);
      writer.Put(var_step, step);
      writer.Put(var_dt, dt);
      writer.Put(var_time, time);
      writer.EndStep();
      writer.Close();
   }

   adios2::IO                    io;
   field_map&                    em_map;
   adios2::Variable<double>      var_ex_energy;
   adios2::Variable<double>      var_ey_energy;
   adios2::Variable<double>      var_ez_energy;
   adios2::Variable<double>      var_bx_energy;
   adios2::Variable<double>      var_by_energy;
   adios2::Variable<double>      var_bz_energy;
   adios2::Variable<std::size_t> var_step;
   adios2::Variable<double>      var_dt;
   adios2::Variable<double>      var_time;
}; // end struct EMTotalEnergyMetric

// =========================================
// ======== Particle Energy Metric =========
struct ParticleTotalEnergyMetric final : detail::MetricBase {
   using group_t = particles::ParticleGroup;
   using group_map = std::unordered_map<std::string, group_t&>;

   struct GroupVariable {
      const group_t& group;
      adios2::Variable<double> variable;
   };

   ParticleTotalEnergyMetric(const std::vector<group_t>& groups_, adios2::IO&& io_)
   : io(io_),
     var_step(io.DefineVariable<std::size_t>("Step")),
     var_dt(io.DefineVariable<double>("dt")),
     var_time(io.DefineVariable<double>("Time"))
   {
      for (auto& g : groups_) {
         groups.push_back(GroupVariable{
            .group = g,
            .variable = io.DefineVariable<double>(g.name, {1}, {0}, {1}, adios2::ConstantDims)
         });
      }

      io.DefineAttribute<std::string>("Unit", "s", "Time");
      io.DefineAttribute<std::string>("Unit", "s", "dt");
   }

   static void write(const auto&, const auto&, const auto, const auto) requires(!(push_enabled and coll_enabled)) {}

   void write(const std::string& dir, const std::string&, const std::size_t step, const double time) override {
      const std::string file{dir + "/particles_energy.bp"};

      if (step % particle_save_interval != 0) { return; }

      adios2::Engine writer = io.Open(file, (step == 0) ? adios2::Mode::Write : adios2::Mode::Append);
      writer.BeginStep();

      for (const auto& [g, var]: groups) {
         auto result = 0.0;
         #pragma omp parallel num_threads(nThreads) default(shared)
         {
            #pragma omp for reduction(+:result)
            for (std::size_t i = 0; i < g.num_particles(); i++) {
               result += (g.particles[i].gamma - 1.0) * g.particles[i].weight;
            }
         }
         result *= g.mass * constants::c_sqr<double>;
         writer.Put(var, result);
      }

      writer.Put(var_step, step);
      writer.Put(var_dt, dt);
      writer.Put(var_time, time);
      writer.EndStep();
      writer.Close();
   }

   adios2::IO                    io;
   std::vector<GroupVariable>    groups;
   adios2::Variable<std::size_t> var_step;
   adios2::Variable<double>      var_dt;
   adios2::Variable<double>      var_time;
}; // end struct ParticleTotalEnergyMetric

// =========================================
// ========= Particle Dump Metric ==========
struct ParticleDumpMetric final : detail::MetricBase {
   using group_t = particles::ParticleGroup;

   ParticleDumpMetric(const group_t& group_, adios2::IO&& io_)
   : io(io_),
     group(group_),
     var_loc(io.DefineVariable<double>("Position", {group.num_particles(), 3}, {0, 0}, {group.num_particles(), 3})),
     var_vel(io.DefineVariable<double>("Velocity", {group.num_particles(), 3}, {0, 0}, {group.num_particles(), 3})),
     var_w(io.DefineVariable<double>("Weight", {group.num_particles(), 1}, {0, 0}, {group.num_particles(), 1})),
     var_gamma(io.DefineVariable<double>("Gamma", {group.num_particles(), 1}, {0, 0}, {group.num_particles(), 1})),
     var_step(io.DefineVariable<std::size_t>("Step", {1}, {0}, {1}, adios2::ConstantDims)),
     var_dt(io.DefineVariable<double>("dt", {1}, {0}, {1}, adios2::ConstantDims)),
     var_time(io.DefineVariable<double>("Time", {1}, {0}, {1}, adios2::ConstantDims))
   {
      io.DefineAttribute<std::string>("Name", group.name);
      io.DefineAttribute<double>("Mass", group.mass);
      io.DefineAttribute<std::string>("Mass/Unit", "kg");
      io.DefineAttribute<double>("Charge", group.charge);
      io.DefineAttribute<std::string>("Charge/Unit", "C");
      io.DefineAttribute<std::size_t>("Atomic Number", group.atomic_number);
      io.DefineAttribute<std::string>("Unit", "m", "Position");
      io.DefineAttribute<std::string>("Unit", "m/s", "Velocity");
      // io.DefineAttribute<std::size_t>("Tracer", group.tracer);
      // io.DefineAttribute<std::size_t>("Sourcer", group.sourcer);
      io.DefineAttribute<std::string>("Unit", "s", "Time");
      io.DefineAttribute<std::string>("Unit", "s", "dt");
      io.DefineAttribute<double>("Cell Volume", dx * dy * dz);
   }

   static void write(const auto&, const auto&, const auto, const auto) requires(!(push_enabled and coll_enabled)) {}

   void write(const std::string& dir, const std::string& step_ext, const std::size_t step, const double time) override {
      static constexpr vec3 delta{dx, dy, dz};
      static constexpr vec3 lb{x_range[0], y_range[0], z_range[0]};

      const std::string file{dir + "/" + group.name + "_dump" + "_" + step_ext};

      if (step % particle_save_interval != 0) { return; }

      const auto nParticles = group.num_particles();
      var_loc.SetShape({nParticles, 3});
      var_loc.SetSelection({{0, 0}, {nParticles, 3}}); // {{start}, {count}}

      var_vel.SetShape({nParticles, 3});
      var_vel.SetSelection({{0, 0}, {nParticles, 3}}); // {{start}, {count}}

      var_w.SetShape({nParticles, 1});
      var_w.SetSelection({{0, 0}, {nParticles, 1}}); // {{start}, {count}}

      var_gamma.SetShape({nParticles, 1});
      var_gamma.SetSelection({{0, 0}, {nParticles, 1}}); // {{start}, {count}}
      
      for (const auto& p: group.particles) {
         for (std::size_t d = 0; d < 3; d++) {
            position.push_back(lb[d] + delta[d] * p.location[d]);
            velocity.push_back(p.velocity[d]);
         }
         weight.push_back(p.weight);
         gamma.push_back(p.gamma);
      }
      
      adios2::Engine writer = io.Open(file, adios2::Mode::Write);
      writer.BeginStep();
      
      writer.Put(var_loc, position.data());
      writer.Put(var_vel, velocity.data());
      writer.Put(var_w, weight.data());
      writer.Put(var_gamma, gamma.data());
      writer.Put(var_step, step);
      writer.Put(var_dt, dt);
      writer.Put(var_time, time);

      writer.EndStep();
      writer.Close();

      position.clear();
      velocity.clear();
      weight.clear();
      gamma.clear();
   }

   adios2::IO                    io;
   const group_t&                group;
   adios2::Variable<double>      var_loc;
   adios2::Variable<double>      var_vel;
   adios2::Variable<double>      var_w;
   adios2::Variable<double>      var_gamma;
   adios2::Variable<std::size_t> var_step;
   adios2::Variable<double>      var_dt;
   adios2::Variable<double>      var_time;
   std::vector<double>           position;
   std::vector<double>           velocity;
   std::vector<double>           weight;
   std::vector<double>           gamma;
}; // end struct ParticleDumpMetric

// // ========================================================
// // ======== Particle Density/Temperature Metric ===========
// struct ParticleMetric final : detail::MetricBase {
//    using group_t = particles::ParticleGroup;
//
//    ParticleMetric(const group_t&    g_,
//                   adios2::IO&&      io_,
//                   const std::size_t ncx,
//                   const std::size_t ncy,
//                   const std::size_t ncz)
//    : io(io_),
//      group(g_),
//      var_density(io.DefineVariable<double>("Density", {ncx, ncy, ncz}, {0, 0, 0}, {ncx, ncy, ncz}, adios2::ConstantDims)),
//      var_temp(io.DefineVariable<double>("Temperature", {ncx, ncy, ncz}, {0, 0, 0}, {ncx, ncy, ncz}, adios2::ConstantDims)),
//      var_step(io.DefineVariable<std::size_t>("Step")),
//      var_dt(io.DefineVariable<double>("dt")),
//      var_time(io.DefineVariable<double>("Time")),
//      density(ncx, ncy, ncz),
//      T_avg(ncx * ncy * ncz),
//      KE_total(ncx * ncy * ncz) {
//       io.DefineAttribute<std::string>("Name", group.name);
//       io.DefineAttribute<double>("Mass", group.mass);
//       io.DefineAttribute<std::string>("Mass/Unit", "kg");
//       io.DefineAttribute<double>("Charge", group.charge);
//       io.DefineAttribute<std::string>("Charge/Unit", "C");
//       io.DefineAttribute<std::size_t>("Atomic Number", group.atomic_number);
//       io.DefineAttribute<double>("Cell Volume", dx * dy * dz);
//       io.DefineAttribute<std::size_t>("dims", dims.data(), 3);
//       io.DefineAttribute<double>("x_range", x_range.data(), 2);
//       io.DefineAttribute<double>("y_range", y_range.data(), 2);
//       io.DefineAttribute<double>("z_range", z_range.data(), 2);
//       io.DefineAttribute<std::string>("File Type", "Metric");
//    }
//
//    void update_metrics() {
//       //using Shape = interp::InterpolationShape<1>::Type;
//       static constexpr auto V_cell_inv = 1.0 / (dx * dy * dz);
//       static constexpr auto temp_coef  = 2.0 / (3.0 * constants::q_e<double>);
//
//       const auto mc2 = group.mass * constants::c_sqr<double>;
//
//       std::ranges::fill(density, 0.0);
//       std::ranges::fill(T_avg, 0.0);
//       std::ranges::fill(KE_total, 0.0);
//
//       // first order density
//       #pragma omp parallel num_threads(nThreads) default(none) shared(mc2)
//       {
//       //    #pragma omp for
//       //    for (std::size_t pid = 0; pid < group->num_particles(); pid++) {
//       //       const auto& p   = group->particles[pid];
//       //       if (p.disabled) { continue; }
//       //       const vec3 loc_half = particles::getCellIndices<float>(p.location + 1.0f) + 0.5f;
//       //       const vec3 hid = loc_half.as_type<std::size_t>();
//       //
//       //       const vec3 p_half = p.location - loc_half;
//       //       const auto shapeI0 = Shape::shape_array(p_half[0]);
//       //       const auto shapeJ0 = Shape::shape_array(p_half[1]);
//       //       const auto shapeK0 = Shape::shape_array(p_half[2]);
//       //
//       //       for (int i = Shape::Begin; i <= Shape::End; ++i) {
//       //          const auto& s0i = shapeI0[i - Shape::Begin];
//       //          for (int j = Shape::Begin; j <= Shape::End; ++j) {
//       //             const auto& s0j = shapeJ0[j - Shape::Begin];
//       //             for (int k = Shape::Begin; k <= Shape::End; ++k) {
//       //                const auto& s0k = shapeK0[k - Shape::Begin];
//       //                #pragma omp atomic update
//       //                density(hid[0] + i, hid[1] + j, hid[2] + k) += s0i * s0j * s0k * p.weight;
//       //             } // end for(k)
//       //          } // end for(j)
//       //       } // end for(i)
//       //    }
//
//
//          #pragma omp for simd
//          for (std::size_t pid = 0; pid < group.num_particles(); pid++) {
//             const auto& p   = group.particles[pid];
//             const auto cid = particles::getCellIndex(p.location);
//             #pragma omp atomic update
//             density[cid] += p.weight;
//             #pragma omp atomic update
//             KE_total[cid] += p.weight * mc2 * (p.gamma - 1.0);
//          }
//
//          #pragma omp for simd
//          for (std::size_t i = 0; i < T_avg.size(); i++) {
//             if (density[i] == 0.0) { continue; }
//             T_avg[i] = temp_coef * KE_total[i] / density[i];
//             //std::cout << KE_total[i] / density[i] << " " << T_avg[i] << std::endl;
//          }
//
//          #pragma omp for simd
//          for (std::size_t i = 0; i < density.size(); i++) {
//             density[i] *= V_cell_inv;
//          }
//       } // end parallel
//    }
//
//    static void write(const auto&, const auto&, const auto, const auto) requires(!particles_enabled) {}
//
//    void write(const std::string& dir, const std::string& step_ext, const std::size_t step, const double time) override {
//       const std::string file{dir + "/" + group.name + "_" + step_ext};
//
//       update_metrics();
//
//       adios2::Engine writer = io.Open(file, adios2::Mode::Write);
//       writer.BeginStep();
//
//       writer.Put(var_density, density.data());
//       writer.Put(var_temp, T_avg.data());
//       writer.Put(var_step, step);
//       writer.Put(var_dt, dt);
//       writer.Put(var_time, time);
//
//       writer.EndStep();
//       writer.Close();
//    } // end write()
//
//    adios2::IO               io;
//    const group_t&           group;
//    adios2::Variable<double> var_density;
//    adios2::Variable<double> var_temp;
//    adios2::Variable<std::size_t> var_step;
//    adios2::Variable<double> var_dt;
//    adios2::Variable<double> var_time;
//    Array3D<double>          density;
//    std::vector<double>      T_avg;
//    std::vector<double>      KE_total;
//    std::array<std::size_t, 3> dims{Nx - 1, Ny - 1, Nz - 1};
// }; // end struct ParticleMetric


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

   void addMetric(std::unique_ptr<detail::MetricBase>&& m) {
      metrics.push_back(std::move(m));
   }

   void add_em_metrics(auto& em) {
      addMetric(
         std::make_unique<EMFieldsMetric>(
            em.emdata.em_map,
            adios.DeclareIO("EMFields")
         )
      );

      // addMetric(
      //    std::make_unique<EMTotalEnergyMetric>(em.emdata.em_map, adios.DeclareIO("EMFieldsTotalEnergy"))
      // );
   }

   void add_particle_metric(const auto& pg) {
      addMetric(
         std::make_unique<ParticleDumpMetric>(pg, adios.DeclareIO("Particles_" + pg.name + (pg.tracer ? "_tracers" : "")))
      );
   }

   void write(const std::size_t step, const double time) const {
      static constexpr size_t padding      = 10;
      std::string             count_padded = std::to_string(step);
      count_padded.insert(count_padded.begin(), padding - count_padded.length(), '0');
      count_padded += file_ext; // default extension is .bp

      for (const auto& m: metrics) {
         m->write(data_dir, count_padded, step, time);
      }
   }

public:
   adios2::ADIOS adios{};

private:
   std::string file_ext{".bp"};
   std::string data_dir{};
   metrics_vec metrics{};
}; // end class Metrics
} // end namespace tf::metrics

#endif //METRICS_HPP
