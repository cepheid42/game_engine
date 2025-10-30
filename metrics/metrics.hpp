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
struct EMFieldSliceMetric final : detail::MetricBase {
   using pointer_t = Array3D<double>*;
   using emdata_t = electromagnetics::EMData;

   struct FieldVariable {
      pointer_t                   field;
      adios2::Variable<double> variable;
      adios2::Variable<double> x_range_var;
      adios2::Variable<double> y_range_var;
      adios2::Variable<double> z_range_var;
   };

   struct RangeData {
      std::vector<double> xrange;
      std::vector<double> yrange;
      std::vector<double> zrange;
   };

   EMFieldSliceMetric(auto& em_map, const std::vector<std::array<std::string,2>>& components_, const std::vector<std::array<std::size_t,6>>& bounds_, adios2::IO&& io_, const bool separate_files_ = false)
   : components(components_), bounds(bounds_), io(io_),
     var_step(io.DefineVariable<std::size_t>("Step")),
     var_dt(io.DefineVariable<double>("dt")),
     var_time(io.DefineVariable<double>("Time")),
     separate_files(separate_files_)
   {
      std::vector<std::string> names;
      std::array<double,3> offsets = {0, 0, 0};
      for (std::size_t i = 0; i < components.size(); i++) {
         auto& [name, component] = components[i];
         auto& bound = bounds[i];

         fields.push_back({
            .field = em_map[component],
            .variable = io.DefineVariable<double>(
               name,
               {1, bound[1] - bound[0], bound[3] - bound[2], bound[5] - bound[4]}, // shape (global)
               {0, 0, 0, 0},                                                       // start (local)
               {1, bound[1] - bound[0], bound[3] - bound[2], bound[5] - bound[4]}, // count (local)
               adios2::ConstantDims
            ),
            .x_range_var = io.DefineVariable<double>(name + "/XRange", {bound[1]-bound[0]}, {0}, {bound[1]-bound[0]}),
            .y_range_var = io.DefineVariable<double>(name + "/YRange", {bound[3]-bound[2]}, {0}, {bound[3]-bound[2]}),
            .z_range_var = io.DefineVariable<double>(name + "/ZRange", {bound[5]-bound[4]}, {0}, {bound[5]-bound[4]}),

         });
         fields_data.emplace_back(bound[1] - bound[0], bound[3] - bound[2], bound[5] - bound[4]);
         std::string component_name = component;
         if (component_name[0] == 'H') { component_name[0] = 'B'; }
         io.DefineAttribute<std::string>("Component", component_name, name);

         if (component[0] == 'H') {
               io.DefineAttribute<std::string>("Unit", "T", name);
               offsets = {component[1] != 'x' ? 0.5 : 0.0, component[1] != 'y' ? 0.5 : 0.0, component[1] != 'z' ? 0.5 : 0.0};
            } else {
               if (component[0] == 'E') { io.DefineAttribute<std::string>("Unit", "V/m", name); }
               else if (component[0] == 'J') { io.DefineAttribute<std::string>("Unit", "A/m^{2}", name); }
               offsets = {component[1] == 'x' ? 0.5 : 0.0, component[1] == 'y' ? 0.5 : 0.0, component[1] == 'z' ? 0.5 : 0.0};
            }

         io.DefineAttribute<std::string>("Unit", "m", name + "/XRange");
         io.DefineAttribute<std::string>("Unit", "m", name + "/YRange");
         io.DefineAttribute<std::string>("Unit", "m", name + "/ZRange");

         names.push_back(name);
         RangeData temp;

         for (std::size_t j = bound[0]; j < bound[1]; j++) {
            temp.xrange.push_back(x_range[0] + (static_cast<double>(j) + offsets[0]) * dx);
         }
         for (std::size_t j = bound[2]; j < bound[3]; j++) {
            temp.yrange.push_back(y_range[0] + (static_cast<double>(j) + offsets[1]) * dy);
         }
         for (std::size_t j = bound[4]; j < bound[5]; j++) {
            temp.zrange.push_back(z_range[0] + (static_cast<double>(j) + offsets[2]) * dz);
         }
         range_datas.push_back(temp);

      }
      io.DefineAttribute<std::string>("Unit", "s", "Time");
      io.DefineAttribute<std::string>("Unit", "s", "dt");
      io.DefineAttribute<std::string>("Cell Volume/Unit", "m^{3}");
      io.DefineAttribute<std::string>("Name", names.data(), names.size());
      io.DefineAttribute<std::string>("File Type", "Field");

   }

   static void write(const auto&, const auto&, const auto, const auto) requires(!fields_enabled) {}

   void write(const std::string& dir, const std::string& step_ext, const std::size_t step, const double time) override {

      const std::string file = separate_files ? dir + "/fields_" + step_ext : dir + "/fields_slices.bp";
      static constexpr auto cell_volume = dx * dy * dz;

      io.DefineAttribute<double>("Cell Volume", cell_volume);

      // append after step 0
      adios2::Engine writer = io.Open(file, (step == 0 or separate_files) ? adios2::Mode::Write : adios2::Mode::Append);
      writer.BeginStep();

      for (std::size_t l = 0; l < bounds.size(); l++) {
         auto& bound = bounds[l];
         auto& [field, variable, x_range_var, y_range_var, z_range_var] = fields[l];
         auto& field_data = fields_data[l];

         if (writer.OpenMode() == adios2::Mode::Write) {
            const auto& [xrange, yrange, zrange] = range_datas[l];
            writer.Put(x_range_var, xrange.data());
            writer.Put(y_range_var, yrange.data());
            writer.Put(z_range_var, zrange.data());
         }

         // update data for this field slice
         for (std::size_t i = 0; i < field_data.nx(); i++) {
            for (std::size_t j = 0; j < field_data.ny(); j++) {
               for (std::size_t k = 0; k < field_data.nz(); k++) {
                  field_data(i, j, k) = field->operator()(i + bound[0], j + bound[2], k + bound[4]) * (components[l][1][0] == 'H' ? constants::mu0<double> : 1.0);
               }
            }
         }

         writer.Put(variable, field_data.data());
         writer.Put(var_step, step);
         writer.Put(var_dt, dt);
         writer.Put(var_time, time);
      }

      writer.EndStep();
      writer.Close();
   }

   std::vector<std::array<std::string,2>> components;
   std::vector<std::array<std::size_t,6>> bounds;
   adios2::IO                    io;
   adios2::Variable<std::size_t> var_step;
   adios2::Variable<double>      var_dt;
   adios2::Variable<double>      var_time;
   std::vector<FieldVariable>    fields;
   std::vector<Array3D<double>>  fields_data;
   std::vector<RangeData>        range_datas;
   bool                          separate_files;
}; // end struct EMFieldSliceMetric

// =====================================
// ======== EM Energy Metric ===========
struct EMTotalEnergyMetric final : detail::MetricBase {
   using pointer_t = Array3D<double>*;
   using field_map = std::unordered_map<std::string, pointer_t>;

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
      io.DefineAttribute<std::string>("Cell Volume/Unit", "m^{3}");
      io.DefineAttribute<std::string>("Unit", "J", "Ex Energy");
      io.DefineAttribute<std::string>("Unit", "J", "Ey Energy");
      io.DefineAttribute<std::string>("Unit", "J", "Ez Energy");
      io.DefineAttribute<std::string>("Unit", "J", "Bx Energy");
      io.DefineAttribute<std::string>("Unit", "J", "By Energy");
      io.DefineAttribute<std::string>("Unit", "J", "Bz Energy");

      const std::vector<std::string> names = {"Ex Energy", "Ey Energy", "Ez Energy", "Bx Energy", "By Energy", "Bz Energy"};
      io.DefineAttribute<std::string>("Name", names.data(), names.size());
      io.DefineAttribute<std::string>("File Type", "Metric");
   }

   static void write(const auto&, const auto&, const auto, const auto) requires(!fields_enabled) {}

   void write(const std::string& dir, const std::string&, const std::size_t step, const double time) override {

      const std::string file{dir + "/fields_energy.bp"};
      static constexpr auto cell_volume = dx * dy * dz;

      io.DefineAttribute<double>("Cell Volume", cell_volume);

      Ex_energy = 0.0;
      Ey_energy = 0.0;
      Ez_energy = 0.0;
      Bx_energy = 0.0;
      By_energy = 0.0;
      Bz_energy = 0.0;

      // todo: fix this to use the right volumes for non-uniform mesh!
#pragma omp parallel num_threads(nThreads) default(none) shared(Ex_energy, Ey_energy, Ez_energy, Bx_energy, By_energy, Bz_energy)
      {
#pragma omp for reduction(+:Ex_energy)
         for (std::size_t i = 0; i<em_map["Ex"]->size(); i++) {
           Ex_energy += math::SQR(em_map["Ex"]->operator[](i)) * cell_volume;
        }
#pragma omp for reduction(+:Ey_energy)
         for (std::size_t i = 0; i<em_map["Ey"]->size(); i++) {
            Ey_energy += math::SQR(em_map["Ey"]->operator[](i)) * cell_volume;
         }
#pragma omp for reduction(+:Ez_energy)
         for (std::size_t i = 0; i<em_map["Ez"]->size(); i++) {
            Ez_energy += math::SQR(em_map["Ez"]->operator[](i)) * cell_volume;
         }
#pragma omp for reduction(+:Bx_energy)
         for (std::size_t i=0; i<em_map["Hx"]->size(); i++) {
            Bx_energy += math::SQR(em_map["Hx"]->operator[](i)) * cell_volume;
         }
#pragma omp for reduction(+:By_energy)
         for (std::size_t i=0; i<em_map["Hy"]->size(); i++) {
            By_energy += math::SQR(em_map["Hy"]->operator[](i)) * cell_volume;
         }
#pragma omp for reduction(+:Bz_energy)
         for (std::size_t i=0; i<em_map["Hz"]->size(); i++) {
            Bz_energy += math::SQR(em_map["Hz"]->operator[](i)) * cell_volume;
         }
      }

      // append after step 0
      adios2::Engine writer = io.Open(file, (step == 0) ? adios2::Mode::Write : adios2::Mode::Append);
      writer.BeginStep();
      writer.Put(var_ex_energy, Ex_energy * tf::constants::one_half<double> * tf::constants::eps0<double>);
      writer.Put(var_ey_energy, Ey_energy * tf::constants::one_half<double> * tf::constants::eps0<double>);
      writer.Put(var_ez_energy, Ez_energy * tf::constants::one_half<double> * tf::constants::eps0<double>);
      writer.Put(var_bx_energy, Bx_energy * tf::constants::one_half<double> * tf::constants::mu0<double>);
      writer.Put(var_by_energy, By_energy * tf::constants::one_half<double> * tf::constants::mu0<double>);
      writer.Put(var_bz_energy, Bz_energy * tf::constants::one_half<double> * tf::constants::mu0<double>);
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
   double Ex_energy = 0.0;
   double Ey_energy = 0.0;
   double Ez_energy = 0.0;
   double Bx_energy = 0.0;
   double By_energy = 0.0;
   double Bz_energy = 0.0;
   std::size_t idx = 0;
}; // end struct EMTotalEnergyMetric

// =========================================
// ======== Particle Energy Metric =========
struct ParticleTotalEnergyMetric final : detail::MetricBase {
   using group_t = particles::ParticleGroup;

   ParticleTotalEnergyMetric(const std::vector<group_t>& groups_, adios2::IO&& io_)
   : io(io_),
     groups(groups_),
     //var_energy(io.DefineVariable<double>("Particle Energy", {groups.size()}, {0}, {groups.size()})),
     var_step(io.DefineVariable<std::size_t>("Step")),
     var_dt(io.DefineVariable<double>("dt")),
     var_time(io.DefineVariable<double>("Time")) {
      std::vector<std::string> names;
      for (auto& g : groups) {
         var_particle_energy.emplace_back(io.DefineVariable<double>(g.name));
         io.DefineAttribute<std::string>("Unit", "J", g.name);
         names.push_back(g.name);
      }
      io.DefineAttribute<std::string>("Unit", "s", "Time");
      io.DefineAttribute<std::string>("Unit", "s", "dt");
      io.DefineAttribute<std::string>("Name", names.data(), names.size());
      particle_energy.reserve(groups.size());

      io.DefineAttribute<std::string>("File Type", "Metric");
   }

   static void write(const auto&, const auto&, const auto, const auto) requires(!particles_enabled) {}

   void write(const std::string& dir, const std::string&, const std::size_t step, const double time) override {

      for (const auto& group: groups) {

         double result = 0.0;

#pragma omp parallel num_threads(nThreads) default(none) shared(group, result)
         {
#pragma omp for reduction(+:result)
            for (std::size_t i=0; i<group.num_particles(); i++) {
               result += (group.particles[i].gamma - 1.0) * group.particles[i].weight;
            }
         }
         particle_energy.push_back(result * group.mass * constants::c_sqr<double>);
      }

      const std::string file{dir + "/particles_energy.bp"};

      adios2::Engine writer = io.Open(file, (step == 0) ? adios2::Mode::Write : adios2::Mode::Append);
      writer.BeginStep();
      for (std::size_t i = 0; i < groups.size(); i++) {
         writer.Put(var_particle_energy[i], particle_energy[i]);
      }
      // writer.Put(var_energy, energy.data());
      writer.Put(var_step, step);
      writer.Put(var_dt, dt);
      writer.Put(var_time, time);
      writer.EndStep();
      writer.Close();

      particle_energy.clear();
   }

   adios2::IO                    io;
   const std::vector<group_t>&   groups;
   adios2::Variable<std::size_t> var_step;
   adios2::Variable<double>      var_dt;
   adios2::Variable<double>      var_time;
   std::vector<adios2::Variable<double>> var_particle_energy;
   // adios2::Variable<double>      var_energy;
   std::vector<double>           particle_energy;
}; // end struct ParticleTotalEnergyMetric

// =========================================
// ========= Particle Dump Metric ==========
struct ParticleDumpMetric final : detail::MetricBase {
   using group_t = particles::ParticleGroup;

   ParticleDumpMetric(const group_t& group_, adios2::IO&& io_, const std::size_t chunk_size_ = 20)
      : io(io_),
        group(group_),
        chunk_size(chunk_size_),
        var_loc(io.DefineVariable<double>("Position", {chunk_size,group.num_particles(), 3}, {0, 0, 0}, {chunk_size,group.num_particles(), 3})),
        var_vel(io.DefineVariable<double>("Velocity", {chunk_size,group.num_particles(), 3}, {0, 0, 0}, {chunk_size,group.num_particles(), 3})),
        var_w(io.DefineVariable<double>("Weight", {chunk_size,group.num_particles()}, {0, 0}, {chunk_size,group.num_particles()})),
        var_gamma(io.DefineVariable<double>("Gamma", {chunk_size,group.num_particles()}, {0, 0}, {chunk_size,group.num_particles()})),
        var_step(io.DefineVariable<std::size_t>("Step", {chunk_size}, {0}, {chunk_size})),
        var_dt(io.DefineVariable<double>("dt", {chunk_size}, {0}, {chunk_size})),
        var_time(io.DefineVariable<double>("Time", {chunk_size}, {0}, {chunk_size}))
   {
      position.reserve(3 * group.num_particles() * chunk_size);
      velocity.reserve(3 * group.num_particles() * chunk_size);
      weight.reserve(group.num_particles() * chunk_size);
      gamma.reserve(group.num_particles() * chunk_size);
      steps.reserve(chunk_size);
      dts.reserve(chunk_size);
      times.reserve(chunk_size);

      io.DefineAttribute<std::string>("Name", group.name);
      io.DefineAttribute<double>("Mass", group.mass);
      io.DefineAttribute<std::string>("Mass/Unit", "kg");
      io.DefineAttribute<double>("Charge", group.charge);
      io.DefineAttribute<std::string>("Charge/Unit", "C");
      io.DefineAttribute<std::size_t>("Atomic Number", group.atomic_number);
      io.DefineAttribute<std::string>("Unit", "m", "Position");
      io.DefineAttribute<std::string>("Unit", "m/s", "Velocity");
      io.DefineAttribute<std::size_t>("Tracer", group.tracer);
      io.DefineAttribute<std::size_t>("Sourcer", group.sourcer);
      io.DefineAttribute<std::string>("Unit", "s", "Time");
      io.DefineAttribute<std::string>("Unit", "s", "dt");
      io.DefineAttribute<std::string>("File Type", "Particle");
   }

   static void write(const auto&, const auto&, const auto, const auto) requires(!particles_enabled) {}
   void write(const std::string& dir, const std::string& step_ext, const std::size_t step, const double time) override {

      static constexpr vec3 delta{dx, dy, dz};
      static constexpr vec3 lb{x_range[0], y_range[0], z_range[0]};

      const std::string file{dir + "/" + group.lower_name + "_dump" + (chunk_size == 1 ? "_" + step_ext : ".bp")};

      // fill tracer data storage
      for (const auto& p: group.particles) {
         for (std::size_t d = 0; d < 3; d++) {
            position.push_back(lb[d] + delta[d] * p.location[d]);
            velocity.push_back(p.velocity[d]);
         }
         weight.push_back(p.weight);
         gamma.push_back(p.gamma);
      }
      steps.push_back(step);
      dts.push_back(dt);
      times.push_back(time);

      idx++;
      if (idx == chunk_size) {
         adios2::Engine writer = io.Open(file, (step == chunk_size - 1 || chunk_size == 1) ? adios2::Mode::Write : adios2::Mode::Append);
         writer.BeginStep();
         writer.Put(var_loc, position.data());
         writer.Put(var_vel, velocity.data());
         writer.Put(var_w, weight.data());
         writer.Put(var_gamma, gamma.data());
         writer.Put(var_step, steps.data());
         writer.Put(var_dt, dts.data());
         writer.Put(var_time, times.data());

         writer.EndStep();
         writer.Close();

         steps.clear();
         dts.clear();
         times.clear();
         position.clear();
         velocity.clear();
         weight.clear();
         gamma.clear();

         idx = 0;
      }
   }

   adios2::IO                    io;
   const group_t&                group;
   std::size_t                   chunk_size;
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
   std::vector<std::size_t>      steps;
   std::vector<double>           dts;
   std::vector<double>           times;
   std::size_t                   idx = 0;
}; // end struct ParticleDumpMetric

// ========================================================
// ======== Particle Density/Temperature Metric ===========
struct ParticleMetric final : detail::MetricBase {
   using group_t = particles::ParticleGroup;

   ParticleMetric(const group_t&    g_,
                  adios2::IO&&      io_,
                  const std::size_t ncx,
                  const std::size_t ncy,
                  const std::size_t ncz)
   : io(io_),
     group(g_),
     var_density(io.DefineVariable<double>("Density", {ncx, ncy, ncz}, {0, 0, 0}, {ncx, ncy, ncz}, adios2::ConstantDims)),
     var_temp(io.DefineVariable<double>("Temperature", {ncx, ncy, ncz}, {0, 0, 0}, {ncx, ncy, ncz}, adios2::ConstantDims)),
     var_step(io.DefineVariable<std::size_t>("Step")),
     var_dt(io.DefineVariable<double>("dt")),
     var_time(io.DefineVariable<double>("Time")),
     density(ncx, ncy, ncz),
     T_avg(ncx * ncy * ncz),
     KE_total(ncx * ncy * ncz) {
      io.DefineAttribute<std::string>("Name", group.name);
      io.DefineAttribute<double>("Mass", group.mass);
      io.DefineAttribute<std::string>("Mass/Unit", "kg");
      io.DefineAttribute<double>("Charge", group.charge);
      io.DefineAttribute<std::string>("Charge/Unit", "C");
      io.DefineAttribute<std::size_t>("Atomic Number", group.atomic_number);

      // todo: adapt this for nonuniform cell volumes
      io.DefineAttribute<double>("Cell Volume", dx * dy * dz);

      io.DefineAttribute<std::size_t>("dims", dims.data(), 3);
      io.DefineAttribute<double>("x_range", x_range.data(), 2);
      io.DefineAttribute<double>("y_range", y_range.data(), 2);
      io.DefineAttribute<double>("z_range", z_range.data(), 2);
      io.DefineAttribute<std::string>("File Type", "Metric");
   }

   void update_metrics() {
      //using Shape = interp::InterpolationShape<1>::Type;
      static constexpr auto V_cell_inv = 1.0 / (dx * dy * dz);
      static constexpr auto temp_coef  = 2.0 / (3.0 * constants::q_e<double>);

      const auto mc2 = group.mass * constants::c_sqr<double>;

      std::ranges::fill(density, 0.0);
      std::ranges::fill(T_avg, 0.0);
      std::ranges::fill(KE_total, 0.0);

      // first order density
      #pragma omp parallel num_threads(nThreads) default(none) shared(mc2)
      {
      //    #pragma omp for
      //    for (std::size_t pid = 0; pid < group->num_particles(); pid++) {
      //       const auto& p   = group->particles[pid];
      //       if (p.disabled) { continue; }
      //       const vec3 loc_half = particles::getCellIndices<float>(p.location + 1.0f) + 0.5f;
      //       const vec3 hid = loc_half.as_type<std::size_t>();
      //
      //       const vec3 p_half = p.location - loc_half;
      //       const auto shapeI0 = Shape::shape_array(p_half[0]);
      //       const auto shapeJ0 = Shape::shape_array(p_half[1]);
      //       const auto shapeK0 = Shape::shape_array(p_half[2]);
      //
      //       for (int i = Shape::Begin; i <= Shape::End; ++i) {
      //          const auto& s0i = shapeI0[i - Shape::Begin];
      //          for (int j = Shape::Begin; j <= Shape::End; ++j) {
      //             const auto& s0j = shapeJ0[j - Shape::Begin];
      //             for (int k = Shape::Begin; k <= Shape::End; ++k) {
      //                const auto& s0k = shapeK0[k - Shape::Begin];
      //                #pragma omp atomic update
      //                density(hid[0] + i, hid[1] + j, hid[2] + k) += s0i * s0j * s0k * p.weight;
      //             } // end for(k)
      //          } // end for(j)
      //       } // end for(i)
      //    }


         #pragma omp for simd
         for (std::size_t pid = 0; pid < group.num_particles(); pid++) {
            const auto& p   = group.particles[pid];
            const auto cid = particles::getCellIndex(p.location);
            #pragma omp atomic update
            density[cid] += p.weight;
            #pragma omp atomic update
            KE_total[cid] += p.weight * mc2 * (p.gamma - 1.0);
         }

         #pragma omp for simd
         for (std::size_t i = 0; i < T_avg.size(); i++) {
            if (density[i] == 0.0) { continue; }
            T_avg[i] = temp_coef * KE_total[i] / density[i];
            //std::cout << KE_total[i] / density[i] << " " << T_avg[i] << std::endl;
         }

         #pragma omp for simd
         for (std::size_t i = 0; i < density.size(); i++) {
            density[i] *= V_cell_inv;
         }
      } // end parallel
   }

   static void write(const auto&, const auto&, const auto, const auto) requires(!particles_enabled) {}

   void write(const std::string& dir, const std::string& step_ext, const std::size_t step, const double time) override {

      const std::string file{dir + "/" + group.lower_name + "_" + step_ext};

      update_metrics();

      adios2::Engine writer = io.Open(file, adios2::Mode::Write);
      writer.BeginStep();

      writer.Put(var_density, density.data());
      writer.Put(var_temp, T_avg.data());
      writer.Put(var_step, step);
      writer.Put(var_dt, dt);
      writer.Put(var_time, time);

      writer.EndStep();
      writer.Close();
   } // end write()

   adios2::IO               io;
   const group_t&           group;
   adios2::Variable<double> var_density;
   adios2::Variable<double> var_temp;
   adios2::Variable<std::size_t> var_step;
   adios2::Variable<double> var_dt;
   adios2::Variable<double> var_time;
   Array3D<double>          density;
   std::vector<double>      T_avg;
   std::vector<double>      KE_total;
   std::array<std::size_t, 3> dims{Nx - 1, Ny - 1, Nz - 1};
}; // end struct ParticleMetric


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

   void addMetric(std::unique_ptr<detail::MetricBase>&& m, const std::size_t interval) { metrics.push_back(std::move(m)); metrics[metrics.size()-1]->interval = interval; }

   void add_em_metrics(auto& em) {
      if constexpr (field_dumps_num > 0) {
         std::vector<std::array<std::string,2>> field_dumps_names_vector;
         std::vector<std::array<std::size_t,6>> field_dumps_bounds_vector;
         for (std::size_t i = 0; i < field_dumps_num; i++) {
            field_dumps_names_vector.push_back(field_dumps_names[i]);
            field_dumps_bounds_vector.push_back(field_dumps_bounds[i]);
         }
         addMetric(
            std::make_unique<EMFieldSliceMetric>(em.emdata.em_map, field_dumps_names_vector, field_dumps_bounds_vector, adios.DeclareIO("EMFieldsDumps"), true),
            save_interval_fields
         );
      }
      addMetric(
         std::make_unique<EMTotalEnergyMetric>(em.emdata.em_map, adios.DeclareIO("EMFieldsTotalEnergy")),
         save_interval_energy
      );
      if constexpr (field_slices_num > 0) {
         std::vector<std::array<std::string,2>> field_slices_names_vector;
         std::vector<std::array<std::size_t,6>> field_slices_bounds_vector;
         for (std::size_t i = 0; i < field_slices_num; i++) {
            field_slices_names_vector.push_back(field_slices_names[i]);
            field_slices_bounds_vector.push_back(field_slices_bounds[i]);
         }
         addMetric(
            std::make_unique<EMFieldSliceMetric>(em.emdata.em_map, field_slices_names_vector, field_slices_bounds_vector, adios.DeclareIO("EMFieldsSlices")),
            save_interval_slices
         );
      }
   }

   void add_particle_metric(const auto& pg) {
      addMetric(
         std::make_unique<ParticleDumpMetric>(pg, adios.DeclareIO("Particles_" + pg.name + (pg.tracer ? "_tracers" : "")), pg.tracer ? 20 : 1),
         pg.tracer ? 1ul : save_interval_particles
            );

      if (!pg.tracer) {
         addMetric(
            std::make_unique<ParticleMetric>(
               pg,
               adios.DeclareIO("Particles_" + pg.name + "_metrics"),
               Nx - 1, Ny - 1, Nz - 1
            ), save_interval_metrics
         );
      }
   }

   void write(const std::size_t step, const double time) const {
      static constexpr size_t padding      = 10;
      std::string             count_padded = std::to_string(step);
      count_padded.insert(count_padded.begin(), padding - count_padded.length(), '0');
      count_padded += file_ext; // default extension is .bp

      for (const auto& m: metrics) {
         if (step % m->interval == 0 or step == Nt) {
            m->write(data_dir, count_padded, step, time);
         }
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
