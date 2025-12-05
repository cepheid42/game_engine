#ifndef METRICS_HPP
#define METRICS_HPP

#include "program_params.hpp"
// #include "array.hpp"
// #include "particles.hpp"
// #include "interpolation.hpp"
// #include "mdspan.hpp"

#include <adios2.h>

#include "particles.hpp"

// #include <unordered_map>
// #include <memory>
// #include <ctime>

namespace tf::metrics {
// =======================================
// ======== Metrics Base Class ===========
// namespace detail {
//    struct MetricBase {
//       virtual ~MetricBase() = default;
//       virtual void write(const std::string&, const std::string&, std::size_t, double) = 0;
//    };
// } // end namespace detail

// template<typename mdspan_t>
// struct FieldMetric {
   // mdspan_t field;
   // adios2::Variable<double> variable;
   // std::vector<double> buffer;

   // void write(auto& writer) const {
   //    std::println("3D Writer {}", writer);
   //    writer.Put(variable, &field.data_handle());
   // }

   // void write(auto& writer) requires (field.rank() == 2) {
   //    std::println("2D Writer");
   //    // for (std::size_t i = 0; i < field.extent(0); ++i) {
   //    //    for (std::size_t j = 0; j < field.extent(1); ++j) {
   //    //       buffer.push_back(field[i, j]);
   //    //    }
   //    // }
   //    // writer.Put(variable, buffer.data());
   //    // buffer.clear();
   // }
   //
   // void write(auto& writer) requires (field.rank() == 1) {
   //    std::println("1D Writer");
   //    // for (std::size_t i = 0; i < field.extent(0); ++i) {
   //    //    buffer.push_back(field[i]);
   //    // }
   //    // writer.Put(variable, buffer.data());
   //    // buffer.clear();
   // }
// };

// =====================================
// ====== EM Field Slice Metric ========
struct EMFieldsMetric {
   static void write(auto& emdata, const std::string& dir, const std::string& step_ext) {
      const std::string file = dir + "/fields_" + step_ext;
      adios2::ADIOS adios{};
      adios2::IO io = adios.DeclareIO("EMFields");

      auto ex_var = io.DefineVariable<double>(
         "Ex",
         {Nx - 1, Ny, Nz}, // shape (global)
         {0, 0, 0},        // start (local)
         {Nx - 1, Ny, Nz}, // count (local)
         adios2::ConstantDims
      );

      auto ey_var = io.DefineVariable<double>(
         "Ey",
         {Nx, Ny - 1, Nz}, // shape (global)
         {0, 0, 0},        // start (local)
         {Nx, Ny - 1, Nz}, // count (local)
         adios2::ConstantDims
      );

      auto ez_var = io.DefineVariable<double>(
         "Ez",
         {Nx, Ny, Nz - 1}, // shape (global)
         {0, 0, 0},        // start (local)
         {Nx, Ny, Nz - 1}, // count (local)
         adios2::ConstantDims
      );

      auto hx_var = io.DefineVariable<double>(
         "Hx",
         {Nx, Ny - 1, Nz - 1}, // shape (global)
         {0, 0, 0},            // start (local)
         {Nx, Ny - 1, Nz - 1}, // count (local)
         adios2::ConstantDims
      );

      auto hy_var = io.DefineVariable<double>(
         "Hy",
         {Nx - 1, Ny, Nz - 1}, // shape (global)
         {0, 0, 0},            // start (local)
         {Nx - 1, Ny, Nz - 1}, // count (local)
         adios2::ConstantDims
      );

      auto hz_var = io.DefineVariable<double>(
         "Hz",
         {Nx - 1, Ny - 1, Nz}, // shape (global)
         {0, 0, 0},            // start (local)
         {Nx - 1, Ny - 1, Nz}, // count (local)
         adios2::ConstantDims
      );

      adios2::Engine writer = io.Open(file, adios2::Mode::Write);
      writer.BeginStep();
      writer.Put(ex_var, emdata.Ex.data());
      writer.Put(ey_var, emdata.Ey.data());
      writer.Put(ez_var, emdata.Ez.data());
      writer.Put(hx_var, emdata.Hx.data());
      writer.Put(hy_var, emdata.Hy.data());
      writer.Put(hz_var, emdata.Hz.data());
      writer.EndStep();
      writer.Close();
   }

}; // end struct EMFieldsMetric

// =========================================
// ========= Particle Dump Metric ==========
struct ParticleDumpMetric  {
   static void write(const auto& grp, const std::string& dir, const std::string& step_ext) {
      static constexpr auto V_cell_inv = 1.0 / (dx * dy * dz);

      const std::string file = dir + "/" + grp.name + "_dump_" + step_ext;
      adios2::ADIOS adios{};
      adios2::IO io = adios.DeclareIO("ParticleDump");

      const auto density_var = io.DefineVariable<double>(
         "Density",
         {Nx - 1, Ny - 1, Nz - 1}, // shape (global)
         {0, 0, 0},                // start (local)
         {Nx - 1, Ny - 1, Nz - 1}, // count (local)
         adios2::ConstantDims
      );

      std::vector<double> density((Nx - 1) * (Ny - 1) * (Nz - 1));

      for (const auto& p : grp.particles) {
         const auto cid = particles::getCellIndex(p.location);
         density[cid] += p.weight;
      }

      for (auto& i : density) {
         i *= V_cell_inv;
      }

      adios2::Engine writer = io.Open(file, adios2::Mode::Write);
      writer.BeginStep();
      writer.Put(density_var, density.data());
      writer.EndStep();
      writer.Close();
   }
   
   // static void write(auto& grp, const std::string& dir, const std::string& step_ext) {
   //    static constexpr vec3 delta{dx, dy, dz};
   //    static constexpr vec3 lb{x_range[0], y_range[0], z_range[0]};
   //    const std::string file = dir + "/" + grp.name + "_dump_" + step_ext;
   //    adios2::ADIOS adios{};
   //    adios2::IO io = adios.DeclareIO("ParticleDump");
   //
   //    const auto num_particles = grp.particles.size();
   //
   //    const auto vel_var = io.DefineVariable<double>(
   //       "Velocity",
   //       {num_particles, 3}, // shape (global)
   //       {0, 0},             // start (local)
   //       {num_particles, 3}, // count (local)
   //       adios2::ConstantDims
   //    );
   //
   //    const auto pos_var = io.DefineVariable<double>(
   //       "Position",
   //       {num_particles, 3}, // shape (global)
   //       {0, 0},             // start (local)
   //       {num_particles, 3}, // count (local)
   //       adios2::ConstantDims
   //    );
   //
   //    const auto weight_var = io.DefineVariable<double>(
   //       "Weight",
   //       {num_particles, 1}, // shape (global)
   //       {0, 0},             // start (local)
   //       {num_particles, 1}, // count (local)
   //       adios2::ConstantDims
   //    );
   //
   //    const auto gamma_var = io.DefineVariable<double>(
   //       "Gamma",
   //       {num_particles, 1}, // shape (global)
   //       {0, 0},             // start (local)
   //       {num_particles, 1}, // count (local)
   //       adios2::ConstantDims
   //    );
   //
   //    std::vector<double> velocities{};
   //    std::vector<double> positions{};
   //    std::vector<double> weights{};
   //    std::vector<double> gamma{};
   //
   //    for (const auto& p : grp.particles) {
   //       for (auto d = 0; d < 3; ++d) {
   //          velocities.push_back(p.velocity[d]);
   //          positions.push_back(lb[d] + delta[d] * p.location[d]);
   //       }
   //       weights.push_back(p.weight);
   //       gamma.push_back(p.gamma);
   //    }
   //
   //    adios2::Engine writer = io.Open(file, adios2::Mode::Write);
   //    writer.BeginStep();
   //    writer.Put(vel_var, velocities.data());
   //    writer.Put(pos_var, positions.data());
   //    writer.Put(weight_var, weights.data());
   //    writer.Put(gamma_var, gamma.data());
   //    writer.EndStep();
   //    writer.Close();
   // }
}; // end struct ParticleDumpMetric


// // =======================================
// // ======== Metrics Superclass ===========
// class Metrics {
// public:
//    explicit Metrics(std::string data_dir_)
//    : data_dir(std::move(data_dir_))
//    {}
//
//    void addMetric(std::unique_ptr<detail::MetricBase>&& m) {}
//
//    void write(const std::size_t step, const double time) const {
//       static constexpr size_t padding      = 10;
//       std::string             count_padded = std::to_string(step);
//       count_padded.insert(count_padded.begin(), padding - count_padded.length(), '0');
//       count_padded += file_ext; // default extension is .bp
//    }
//
// public:
//    adios2::ADIOS adios{};
//
// private:
//    std::string file_ext{".bp"};
//    std::string data_dir{};
// }; // end class Metrics
} // end namespace tf::metrics

#endif //METRICS_HPP
