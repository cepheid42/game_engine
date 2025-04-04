#include "metrics.hpp"

#include "program_params.hpp"
#include "particles.hpp"

#include <vector>
#include <print>
#include <algorithm>
#include <adios2.h>

namespace tf::metrics {

  EMFieldsMetric::EMFieldsMetric(const field_map& fields_, adios2::IO&& io_)
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

  void EMFieldsMetric::write(const std::string& dir, const std::string& step_ext) {
    const std::string file{dir + "/fields_" + step_ext};
    adios2::Engine writer = io.Open(file, adios2::Mode::Write);
    writer.BeginStep();

    // ReSharper disable once CppUseElementsView
    for (auto& [field, variable]: fields) {
      writer.Put(variable, field->data());
    }

    writer.EndStep();
    writer.Close();
  }

  ParticleDumpMetric::ParticleDumpMetric(const group_t* group_, adios2::IO&& io_)
  : io(io_),
    group(group_),
    var_loc(io.DefineVariable<compute_t>("Position", {group->num_particles, 3}, {0, 0}, {group->num_particles, 3})),
    var_vel(io.DefineVariable<double>("Velocity", {group->num_particles, 3}, {0, 0}, {group->num_particles, 3})),
    var_w(io.DefineVariable<compute_t>("Weight", {group->num_particles, 1}, {0, 0}, {group->num_particles, 1}))
  {}

  void ParticleDumpMetric::write(const std::string& dir, const std::string& step_ext) {
    const std::string file{dir + "/" + group->name + "_dump_" + step_ext};

    // io.DefineAttribute<std::string>("name", group->name);
    // io.DefineAttribute<std::size_t>("num_particles", group->num_particles);
    // io.DefineAttribute<std::size_t>("atomic_number", group->atomic_number);
    // io.DefineAttribute<compute_t>("mass", group->mass);
    // io.DefineAttribute<compute_t>("charge", group->charge);

    adios2::Engine writer = io.Open(file, adios2::Mode::Write);
    writer.BeginStep();

    constexpr vec3 delta{dx, dy, dz};
    constexpr vec3 lb{x_range[0], y_range[0], z_range[0]};
    const auto& nParticles = group->num_particles;
    std::vector<compute_t> position{};
    std::vector<double> velocity{};
    std::vector<compute_t> weight{};

    position.reserve(3 * nParticles);
    velocity.reserve(3 * nParticles);
    weight.reserve(nParticles);

    for (const auto& cell: group->cells) {
      const auto& idx = particles::morton_decode(cell.cid);
      for (const auto& p: cell.particles) {
        for (std::size_t d = 0; d < 3; d++) {
          position.push_back(lb[d] + delta[d] * (static_cast<compute_t>(idx[d]) + p.location[d]));
          velocity.push_back(p.velocity[d]);
        }
        weight.push_back(p.weight);
      }
    }

    writer.Put(var_loc, position.data());
    writer.Put(var_vel, velocity.data());
    writer.Put(var_w, weight.data());

    writer.EndStep();
    writer.Close();
  }

  ParticleMetric::ParticleMetric(const group_t* g_, adios2::IO&& io_)
  : io(io_),
    group(g_),
    var_density(io.DefineVariable<compute_t>("Density", {Nx - 1, Ny - 1, Nz - 1}, {0, 0, 0}, {Nx - 1, Ny - 1, Nz - 1}, adios2::ConstantDims)),
    var_temp(io.DefineVariable<compute_t>("Temperature", {Nx - 1, Ny - 1, Nz - 1}, {0, 0, 0}, {Nx - 1, Ny - 1, Nz - 1}, adios2::ConstantDims)),
    density(g_->cells.size()),
    T_avg(g_->cells.size())
  {}

  void ParticleMetric::update_metrics() {
    // density = cell_weight / cell_volume
    // Te = dv2sum * (me / (qe * cell_weight)) / 3
    constexpr auto V_cell_inv = (1.0_fp / dx) + (1.0_fp / dy) + (1.0_fp / dz);

    for (std::size_t i = 0; i < group->cells.size(); i++) {
      const auto& [particles, cid] = group->cells[i];

      const auto w_cell = std::ranges::fold_left(particles.begin(), particles.end(), 0.0_fp,
        [](const compute_t init, const particles::Particle& p) { return init + p.weight; }
      );

      auto v_avg = std::ranges::fold_left(particles.begin(), particles.end(), vec3{0.0_fp, 0.0_fp, 0.0_fp},
        [](const vec3<compute_t> init, const particles::Particle& p) { return init + p.weight * (p.velocity / static_cast<compute_t>(p.gamma)); }
      );

      v_avg /= w_cell;

      const auto dv2_sum = std::ranges::fold_left(particles.begin(), particles.end(), 0.0_fp,
        [&v_avg](const compute_t init, const particles::Particle& p) { return init + p.weight * ((p.velocity / static_cast<compute_t>(p.gamma)) - v_avg).length_squared(); }
      );

      // for electrons only
      T_avg[i] = dv2_sum * group->mass / (3.0_fp * w_cell * static_cast<compute_t>(constants::q_e));
      density[i] = w_cell * V_cell_inv;
    } // end for(i)
  }

  void ParticleMetric::write([[maybe_unused]] const std::string& dir,[[maybe_unused]]  const std::string& step_ext) {
    update_metrics();

    const std::string file{dir + "/" + group->name + "_" + step_ext};
    adios2::Engine writer = io.Open(file, adios2::Mode::Write);
    writer.BeginStep();

    writer.Put(var_density, density.data());
    writer.Put(var_temp, T_avg.data());

    writer.EndStep();
    writer.Close();
  } // end write()

  void Metrics::write(const std::size_t step) {
    static constexpr size_t padding = 10;
    std::string count_padded = std::to_string(step);
    count_padded.insert(count_padded.begin(), padding - count_padded.length(), '0');
    count_padded += file_ext; // default extension is .bp

    for (const auto& m: metrics) {
      m->write(data_dir, count_padded);
    }
  }
}

