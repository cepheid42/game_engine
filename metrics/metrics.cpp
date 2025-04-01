#include "metrics.hpp"

#include "program_params.hpp"
#include "particles.hpp"

#include <vector>
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

  ParticleMetric::ParticleMetric(const group_t* group_, adios2::IO&& io_)
  : io(io_),
    group(group_),
    var_loc(io.DefineVariable<compute_t>("Position", {group->num_particles, 3}, {0, 0}, {group->num_particles, 3})),
    var_vel(io.DefineVariable<compute_t>("Velocity", {group->num_particles, 3}, {0, 0}, {group->num_particles, 3})),
    var_w(io.DefineVariable<compute_t>("Weight", {group->num_particles, 1}, {0, 0}, {group->num_particles, 1}))
  {}

  void ParticleMetric::write(const std::string& dir, const std::string& step_ext) {
    const std::string file{dir + "/" + group->name + "_" + step_ext};

    io.DefineAttribute<std::string>("name", group->name);
    io.DefineAttribute<std::size_t>("num_particles", group->num_particles);
    io.DefineAttribute<std::size_t>("atomic_number", group->atomic_number);
    io.DefineAttribute<compute_t>("mass", group->mass);
    io.DefineAttribute<compute_t>("charge", group->charge);

    adios2::Engine writer = io.Open(file, adios2::Mode::Write);
    writer.BeginStep();

    const auto& nParticles = group->num_particles;
    std::vector<compute_t> position(nParticles);
    std::vector<compute_t> velocity(nParticles);
    std::vector<compute_t> weight(nParticles);

    for (const auto& cell: group->cells) {
      for (const auto& p: cell.particles) {
        for (std::size_t i = 0; i < 3; i++) {
          position.push_back(p.location[i]);
          velocity.push_back(p.velocity[i]);
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

