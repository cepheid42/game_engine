#ifndef METRICS_HPP
#define METRICS_HPP

#include "array.hpp"
#include "em_data.hpp"
#include "em_definitions.hpp"
#include "interpolation.hpp"
#include "mdspan.hpp"
#include "particles.hpp"
#include "program_params.hpp"

#include <adios2.h>
#include <algorithm>
#include <memory>
#include <print>
#include <unordered_map>

namespace tf::metrics {
// =======================================
// ======== Metrics Base Class ===========
namespace detail {
   struct MetricBase {
      virtual ~MetricBase() = default;
      virtual void write(const std::string&, const std::string&, std::size_t, double) = 0;
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
   {
      for (auto& [name, field] : em_map) {
         auto nm = name;
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
}; // end struct EMFieldsMetric

// =====================================
// ======== EM Energy Metric ===========
struct EMTotalEnergyMetric final : detail::MetricBase {
   using field_map = std::unordered_map<std::string, Array3D<double>&>;

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
      using mdspan_t = std::mdspan<double, std::dextents<std::size_t, 3>, std::layout_stride>;

      static constexpr auto cell_volume = dx * dy * dz;
      const std::string file{dir + "/fields_energy.bp"};

      if (step % em_save_interval != 0) { return; }

      auto Ex_energy = 0.0;
      auto Ey_energy = 0.0;
      auto Ez_energy = 0.0;
      auto Hx_energy = 0.0;
      auto Hy_energy = 0.0;
      auto Hz_energy = 0.0;

      constexpr auto x0 = x_collapsed ? 0 : PMLDepth;
      constexpr auto y0 = y_collapsed ? 0 : PMLDepth;
      constexpr auto z0 = z_collapsed ? 0 : PMLDepth;
      constexpr auto nx = x_collapsed ? 1 : Nx - PMLDepth;
      constexpr auto ny = y_collapsed ? 1 : Ny - PMLDepth;
      constexpr auto nz = z_collapsed ? 1 : Nz - PMLDepth;
      constexpr auto ncx = x_collapsed ? 1 : Nx - PMLDepth - 1;
      constexpr auto ncy = y_collapsed ? 1 : Ny - PMLDepth - 1;
      constexpr auto ncz = z_collapsed ? 1 : Nz - PMLDepth - 1;
      
      const auto Ex = mdspan_t{&em_map.at("Ex")(x0, y0, z0), {std::extents{ncx, ny, nz}, electromagnetics::ex_stride}};
      const auto Ey = mdspan_t{&em_map.at("Ey")(x0, y0, z0), {std::extents{nx, ncy, nz}, electromagnetics::ey_stride}};
      const auto Ez = mdspan_t{&em_map.at("Ez")(x0, y0, z0), {std::extents{nx, ny, ncz}, electromagnetics::ez_stride}};
      const auto Hx = mdspan_t{&em_map.at("Hx")(x0, y0, z0), {std::extents{nx, ncy, ncz}, electromagnetics::hx_stride}};
      const auto Hy = mdspan_t{&em_map.at("Hy")(x0, y0, z0), {std::extents{ncx, ny, ncz}, electromagnetics::hy_stride}};
      const auto Hz = mdspan_t{&em_map.at("Hz")(x0, y0, z0), {std::extents{ncx, ncy, nz}, electromagnetics::hz_stride}};

      #pragma omp parallel num_threads(nThreads) default(shared)
      {
         #pragma omp for collapse(3) reduction(+:Ex_energy)
         for (auto i = 0zu; i < Ex.extent(0); i++) {
            for (auto j = 0zu; j < Ex.extent(1); j++) {
               for (auto k = 0zu; k < Ex.extent(2); k++) {
                  Ex_energy += math::SQR(Ex[i, j, k]);
               }
            }
         }
         #pragma omp for collapse(3) reduction(+:Ey_energy)
         for (auto i = 0zu; i < Ey.extent(0); i++) {
            for (auto j = 0zu; j < Ey.extent(1); j++) {
               for (auto k = 0zu; k < Ey.extent(2); k++) {
                  Ey_energy += math::SQR(Ey[i, j, k]);
               }
            }
         }
         #pragma omp for collapse(3) reduction(+:Ez_energy)
         for (auto i = 0zu; i < Ez.extent(0); i++) {
            for (auto j = 0zu; j < Ez.extent(1); j++) {
               for (auto k = 0zu; k < Ez.extent(2); k++) {
                  Ez_energy += math::SQR(Ez[i, j, k]);
               }
            }
         }

         #pragma omp for collapse(3) reduction(+:Hx_energy)
         for (auto i = 0zu; i < Hx.extent(0); i++) {
            for (auto j = 0zu; j < Hx.extent(1); j++) {
               for (auto k = 0zu; k < Hx.extent(2); k++) {
                  Hx_energy += math::SQR(Hx[i, j, k]);
               }
            }
         }
         #pragma omp for collapse(3) reduction(+:Hy_energy)
         for (auto i = 0zu; i < Hy.extent(0); i++) {
            for (auto j = 0zu; j < Hy.extent(1); j++) {
               for (auto k = 0zu; k < Hy.extent(2); k++) {
                  Hy_energy += math::SQR(Hy[i, j, k]);
               }
            }
         }
         #pragma omp for collapse(3) reduction(+:Hz_energy)
         for (auto i = 0zu; i < Hz.extent(0); i++) {
            for (auto j = 0zu; j < Hz.extent(1); j++) {
               for (auto k = 0zu; k < Hz.extent(2); k++) {
                  Hz_energy += math::SQR(Hz[i, j, k]);
               }
            }
         }
      } // end parallel

      // append after step 0
      adios2::Engine writer = io.Open(file, (step == 0) ? adios2::Mode::Write : adios2::Mode::Append);
      writer.BeginStep();
      writer.Put(var_ex_energy, Ex_energy * 0.5 * constants::eps0 * cell_volume);
      writer.Put(var_ey_energy, Ey_energy * 0.5 * constants::eps0 * cell_volume);
      writer.Put(var_ez_energy, Ez_energy * 0.5 * constants::eps0 * cell_volume);
      writer.Put(var_bx_energy, Hx_energy * 0.5 * constants::mu0 * cell_volume);
      writer.Put(var_by_energy, Hy_energy * 0.5 * constants::mu0 * cell_volume);
      writer.Put(var_bz_energy, Hz_energy * 0.5 * constants::mu0 * cell_volume);
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
   using group_map = std::unordered_map<std::string, group_t>;

   struct GroupVariable {
      const group_t& group;
      adios2::Variable<double> variable;
   };

   ParticleTotalEnergyMetric(const auto& group_map, adios2::IO&& io_)
   : io(io_),
     var_step(io.DefineVariable<std::size_t>("Step")),
     var_dt(io.DefineVariable<double>("dt")),
     var_time(io.DefineVariable<double>("Time"))
   {
      for (auto& [name, g] : group_map) {
         groups.push_back(GroupVariable{
            .group = g,
            .variable = io.DefineVariable<double>(name, {1}, {0}, {1}, adios2::ConstantDims)
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

      for (const auto& gv: groups) {
         const auto& group = gv.group;
         auto result = 0.0;

         #pragma omp parallel for num_threads(nThreads) default(shared) reduction(+:result)
         for (std::size_t i = 0; i < group.num_particles(); i++) {
            const auto& p = group.particles[i];

            const auto [inew, jnew, knew] = particles::getCellIndices(p.location);
            const auto disabled = (!x_collapsed and (inew < PMLDepth or inew > Nx - 2 - PMLDepth)) or
                                  (!y_collapsed and (jnew < PMLDepth or jnew > Ny - 2 - PMLDepth)) or
                                  (!z_collapsed and (knew < PMLDepth or knew > Nz - 2 - PMLDepth));
            
            if (disabled or p.is_disabled()) { continue; }
            result += (p.gamma - 1.0) * static_cast<double>(p.weight);
         }
         result *= group.mass * constants::c_sqr;
         writer.Put(gv.variable, result);
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
     var_gamma(io.DefineVariable<double>("Gamma", {group.num_particles(), 1}, {0, 0}, {group.num_particles(), 1})),
     var_w(io.DefineVariable<float>("Weight", {group.num_particles(), 1}, {0, 0}, {group.num_particles(), 1})),
     var_step(io.DefineVariable<std::size_t>("Step", {1}, {0}, {1}, adios2::ConstantDims)),
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
      io.DefineAttribute<std::size_t>("Tracer", group.is_tracer);
      io.DefineAttribute<std::string>("Unit", "s", "Time");
      io.DefineAttribute<double>("dt", dt);
   }

   static void write(const auto&, const auto&, const auto, const auto) requires(!(push_enabled and coll_enabled)) {}

   void write(const std::string& dir, const std::string& step_ext, const std::size_t step, const double time) override {
      static constexpr vec3 delta{dx, dy, dz};
      static constexpr vec3 lb{x_range[0], y_range[0], z_range[0]};

      const std::string file{dir + "/" + group.name + "_dump_" + step_ext};

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
         if (p.is_disabled()) { continue; }

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
   adios2::Variable<double>      var_gamma;
   adios2::Variable<float>       var_w;
   adios2::Variable<std::size_t> var_step;
   adios2::Variable<double>      var_time;
   std::vector<double>           position{};
   std::vector<double>           velocity{};
   std::vector<double>           gamma{};
   std::vector<float>            weight{};
}; // end struct ParticleDumpMetric

// =========================================
// ========= Particle Dump Metric ==========
struct ParticleTracerMetric final : detail::MetricBase {
   using group_t = particles::ParticleGroup;

   ParticleTracerMetric(const group_t& group_, adios2::IO&& io_)
   : io(io_),
     group(group_),
     var_p(io.DefineVariable<double>("Position", {particle_save_interval, 3}, {0, 0}, {particle_save_interval, 3}, adios2::ConstantDims)),
     var_v(io.DefineVariable<double>("Velocity", {particle_save_interval, 3}, {0, 0}, {particle_save_interval, 3}, adios2::ConstantDims)),
     var_w(io.DefineVariable<double>("Weight",   {particle_save_interval, 1}, {0, 0}, {particle_save_interval, 1}, adios2::ConstantDims)),
     var_g(io.DefineVariable<double>("Gamma",    {particle_save_interval, 1}, {0, 0}, {particle_save_interval, 1}, adios2::ConstantDims)),
     var_t(io.DefineVariable<double>("Time",     {particle_save_interval, 1}, {0, 0}, {particle_save_interval, 1}, adios2::ConstantDims)),
     var_s(io.DefineVariable<std::size_t>("Step",{particle_save_interval, 1}, {0, 0}, {particle_save_interval, 1}, adios2::ConstantDims))
   {
      io.DefineAttribute<double>("Mass", group.mass);
      io.DefineAttribute<double>("Charge", group.charge);
      io.DefineAttribute<double>("dt", dt);
      io.DefineAttribute<std::size_t>("Atomic Number", group.atomic_number);
      io.DefineAttribute<std::size_t>("Tracer", group.is_tracer);
      io.DefineAttribute<std::string>("Name", group.name);
      io.DefineAttribute<std::string>("Mass/Unit", "kg");
      io.DefineAttribute<std::string>("Charge/Unit", "C");
      io.DefineAttribute<std::string>("Unit", "m", "Position");
      io.DefineAttribute<std::string>("Unit", "m/s", "Velocity");
      io.DefineAttribute<std::string>("Unit", "s", "Time");
   }

   void write(const std::string& dir, const std::string& step_ext, const std::size_t step, const double time) override {
      static constexpr std::array delta{dx, dy, dz};
      static constexpr std::array lb{x_range[0], y_range[0], z_range[0]};

      for (const auto& p: group.particles) {
         // if (p.is_disabled()) { continue; }

         for (std::size_t d = 0; d < 3; d++) {
            position.push_back(lb[d] + delta[d] * p.location[d]);
            velocity.push_back(p.velocity[d]);
         }
         weight.push_back(p.weight);
         gamma.push_back(p.gamma);
      }
      times.push_back(time);
      steps.push_back(step);


      if (step % particle_save_interval == 0zu and step > 0zu) {
         const std::string file{dir + "/" + group.name + "_tracer_" + step_ext};
         adios2::Engine writer = io.Open(file, adios2::Mode::Write);

         writer.BeginStep();

         writer.Put(var_p, position.data());
         writer.Put(var_v, velocity.data());
         writer.Put(var_w, weight.data());
         writer.Put(var_g, gamma.data());
         writer.Put(var_t, times.data());
         writer.Put(var_s, steps.data());

         writer.EndStep();
         writer.Close();

         position.clear();
         velocity.clear();
         weight.clear();
         gamma.clear();
         times.clear();
         steps.clear();
      }
   } // end write()

   adios2::IO                    io;
   const group_t&                group;
   adios2::Variable<double>      var_p;
   adios2::Variable<double>      var_v;
   adios2::Variable<double>      var_w;
   adios2::Variable<double>      var_g;
   adios2::Variable<double>      var_t;
   adios2::Variable<std::size_t> var_s;
   std::vector<double>           position{};
   std::vector<double>           velocity{};
   std::vector<double>           weight{};
   std::vector<double>           gamma{};
   std::vector<double>           times{};
   std::vector<std::size_t>      steps{};
}; // end struct ParticleTracerMetric

// ========================================================
// ======== Particle Density/Temperature Metric ===========
struct ParticleMetric final : detail::MetricBase {
   using group_t = particles::ParticleGroup;
   using XShape = interp::InterpolationShape<x_collapsed ? 1 : interpolation_order>::Type;
   using YShape = interp::InterpolationShape<y_collapsed ? 1 : interpolation_order>::Type;
   using ZShape = interp::InterpolationShape<z_collapsed ? 1 : interpolation_order>::Type;

   ParticleMetric(const group_t& g_, adios2::IO&& io_)
   : io(io_),
     group(g_),
     var_density(io.DefineVariable<double>("Density", {Nx - 1, Ny - 1, Nz - 1}, {0, 0, 0},  {Nx - 1, Ny - 1, Nz - 1}, adios2::ConstantDims)),
     var_temp(io.DefineVariable<double>("Temperature",  {Nx - 1, Ny - 1, Nz - 1}, {0, 0, 0},  {Nx - 1, Ny - 1, Nz - 1}, adios2::ConstantDims)),
     var_step(io.DefineVariable<std::size_t>("Step")),
     var_dt(io.DefineVariable<double>("dt")),
     var_time(io.DefineVariable<double>("Time")),
     density(Nx - 1, Ny - 1, Nz - 1),
     T_avg(Nx - 1, Ny - 1,  Nz - 1),
     KE_total(Nx - 1, Ny - 1,  Nz - 1)
   {
     io.DefineAttribute<std::string>("Name", group.name);
     io.DefineAttribute<double>("Mass", group.mass);
     io.DefineAttribute<std::string>("Mass/Unit", "kg");
     io.DefineAttribute<double>("Charge", group.charge);
     io.DefineAttribute<std::string>("Charge/Unit", "C");
     io.DefineAttribute<std::size_t>("Atomic Number", group.atomic_number);
     io.DefineAttribute<double>("Cell Volume", dx * dy * dz);
     io.DefineAttribute<std::size_t>("dims", dims.data(), 3);
     io.DefineAttribute<double>("x_range", x_range.data(), 2);
     io.DefineAttribute<double>("y_range", y_range.data(), 2);
     io.DefineAttribute<double>("z_range", z_range.data(), 2);
     io.DefineAttribute<std::string>("File Type", "Metric");
   }

   void update_metrics() {
      static constexpr auto V_cell_inv = 1.0 / (dx * dy * dz);
      static constexpr auto temp_coef  = 2.0 / (3.0 * constants::q_e);
      static constexpr vec3 offset{
         interpolation_order == 2 ? 0.5 : 0.0,
         interpolation_order == 2 ? 0.5 : 0.0,
         interpolation_order == 2 ? 0.5 : 0.0
      };

      const auto mc2 = group.mass * constants::c_sqr;

      std::ranges::fill(density, 0.0);
      std::ranges::fill(T_avg, 0.0);
      // std::ranges::fill(KE_total, 0.0);

      #pragma omp parallel num_threads(nThreads) default(shared)
      {
         // #pragma omp for
         // for (auto pid = 0zu; pid < group.num_particles(); pid++) {
         //    const auto& p   = group.particles[pid];
         //    if (p.is_disabled()) { continue; }
         //
         //    const vec3 loc_half = particles::getCellIndices<double>(p.location - offset);
         //    const vec3 hid = loc_half.to_uint();
         //
         //    const vec3 p_half = p.location - loc_half;
         //
         //    const auto shapeI0 = XShape::shape_array(p_half[0]);
         //    const auto shapeJ0 = YShape::shape_array(p_half[1]);
         //    const auto shapeK0 = ZShape::shape_array(p_half[2]);
         //
         //    for (auto i = XShape::Begin; i <= XShape::End; ++i) {
         //       const auto& s0i = shapeI0[i];
         //       for (auto j = YShape::Begin; j <= YShape::End; ++j) {
         //          const auto& s0j = shapeJ0[j];
         //          for (auto k = ZShape::Begin; k <= ZShape::End; ++k) {
         //             const auto& s0k = shapeK0[k];
         //             #pragma omp atomic update
         //             density(hid[0] + i, hid[1] + j, hid[2] + k) += s0i * s0j * s0k * p.weight;
         //             #pragma omp atomic update
         //             T_avg(hid[0] + i, hid[1] + j, hid[2] + k) += s0i * s0j * s0k * p.weight * mc2 * (p.gamma - 1.0);
         //          } // end for(k)
         //       } // end for(j)
         //    } // end for(i)
         // }

         #pragma omp for
         for (std::size_t pid = 0; pid < group.num_particles(); pid++) {
            const auto& p   = group.particles[pid];
            const auto [i, j, k] = particles::getCellIndices(p.location);
            #pragma omp atomic update
            density(i, j, k) += p.weight;
            #pragma omp atomic update
            T_avg(i, j, k) += p.weight * mc2 * (p.gamma - 1.0);
         }

         #pragma omp for
         for (auto i = 0zu; i < T_avg.size(); i++) {
            if (density[i] == 0.0) { continue; }
            T_avg[i] *= temp_coef / density[i];
         }

         #pragma omp for
         for (auto i = 0zu; i < density.size(); i++) {
            density[i] *= V_cell_inv;
         }
      } // end parallel
   }

   void write(const std::string& dir, const std::string& step_ext, const std::size_t step, const double time) override {
      const std::string file{dir + "/" + group.name + "_" + step_ext};

      if (step % particle_save_interval != 0) { return; }

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
   Array3D<double>          T_avg;
   Array3D<double>          KE_total;
   std::array<std::size_t, 3> dims{Nx - 1, Ny - 1, Nz - 1};
}; // end struct ParticleMetric

// =======================================
// ======== Metrics Superclass ===========
class Metrics {
   using metrics_vec = std::vector<std::unique_ptr<detail::MetricBase>>;

public:
   explicit Metrics(std::string data_dir_, const auto& spec, auto& em_map, auto& p_groups)
   : data_dir(std::move(data_dir_))
   {
      // const std::time_t time = std::time({});
      // char timestr[std::size("YYMMDDHHMM")];
      // std::strftime(std::data(timestr), std::size(timestr), "%y%m%d%H%M", std::localtime(&time));
      //
      // data_dir += "_" + std::string(timestr);

      for (const auto& s : spec) {
         switch (s) {
            case MetricType::ParticleDump:
               {
                  for (auto& [name, g] : p_groups) {
                     if (g.is_tracer) {
                        metrics.push_back(std::make_unique<ParticleTracerMetric>(g, adios.DeclareIO("Particles_" + name + "_tracer")));
                     } else {
                        metrics.push_back(std::make_unique<ParticleDumpMetric>(g, adios.DeclareIO("Particles_" + name + "_dump")));
                     }
                  }
                  break;
               }
            case MetricType::ParticleDiag:
               {
                  for (auto& [name, g] : p_groups) {
                     metrics.push_back(std::make_unique<ParticleMetric>(g, adios.DeclareIO("Particles_" + name)));
                  }
                  break;
               }
            case MetricType::ParticleEnergy:
               {
                  metrics.push_back(std::make_unique<ParticleTotalEnergyMetric>(p_groups, adios.DeclareIO("Particle_energy")));
                  break;
               }
            case MetricType::FieldDump:
               {
                  metrics.push_back(std::make_unique<EMFieldsMetric>(em_map,  adios.DeclareIO("EMFields_dump")));
                  break;
               }
            case MetricType::FieldEnergy:
               {
                  metrics.push_back(std::make_unique<EMTotalEnergyMetric>(em_map,  adios.DeclareIO("EMFields_energy")));
                  break;
               }
            default:
               assert(false);
               break;
         } // end switch(s)
      } // end for(s)
   } // end Metric ctor


   void write(const std::size_t step, const double time) const {
      static constexpr auto padding = 10zu;
      auto count_padded = std::to_string(step);
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
