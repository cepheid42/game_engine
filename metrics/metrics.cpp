#include "metrics.hpp"

#include "program_params.hpp"
#include "particles.hpp"
#include "constants.hpp"

#include <vector>
// #include <print>
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
    var_vel(io.DefineVariable<compute_t>("Velocity", {group->num_particles, 3}, {0, 0}, {group->num_particles, 3})),
    var_w(io.DefineVariable<compute_t>("Weight", {group->num_particles, 1}, {0, 0}, {group->num_particles, 1})),
    var_gamma(io.DefineVariable<double>("Gamma", {group->num_particles, 1}, {0, 0}, {group->num_particles, 1}))
  {}

  void ParticleDumpMetric::write(const std::string& dir, const std::string& step_ext) {
    const std::string file{dir + "/" + group->name + "_dump_" + step_ext};

    // io.DefineAttribute<std::string>("name", group->name);
    // io.DefineAttribute<std::size_t>("num_particles", group->num_particles);
    // io.DefineAttribute<std::size_t>("atomic_number", group->atomic_number);
    // io.DefineAttribute<compute_t>("mass", group->mass);
    // io.DefineAttribute<compute_t>("charge", group->charge);

    var_loc.SetShape({group->num_particles, 3});
    var_loc.SetSelection({{0, 0}, {group->num_particles, 3}}); // {{start}, {count}}

    var_vel.SetShape({group->num_particles, 3});
    var_vel.SetSelection({{0, 0}, {group->num_particles, 3}}); // {{start}, {count}}

    var_w.SetShape({group->num_particles, 1});
    var_w.SetSelection({{0, 0}, {group->num_particles, 1}}); // {{start}, {count}}

    var_gamma.SetShape({group->num_particles, 1});
    var_gamma.SetSelection({{0, 0}, {group->num_particles, 1}}); // {{start}, {count}}

    adios2::Engine writer = io.Open(file, adios2::Mode::Write);
    writer.BeginStep();

    constexpr vec3 delta{dx, dy, dz};
    constexpr vec3 lb{x_range[0], y_range[0], z_range[0]};
    const auto& nParticles = group->num_particles;

    std::vector<compute_t> position{};
    std::vector<compute_t> velocity{};
    std::vector<compute_t> weight{};
    std::vector<double> gamma{};

    position.reserve(3 * nParticles);
    velocity.reserve(3 * nParticles);
    weight.reserve(nParticles);
    gamma.reserve(nParticles);

    for (const auto& [chunks, idxs]: group->cells) {
      for (const auto& chunk : chunks) {
        for (std::size_t pid = 0; pid < particles::ParticleChunk::n_particles; pid++) {
          if (!chunk.active.test(pid)) { continue; }
          const auto& p = chunk[pid];
          for (std::size_t d = 0; d < 3; d++) {
            position.push_back(lb[d] + delta[d] * (static_cast<compute_t>(idxs[d]) + p.location[d]));
            velocity.push_back(p.velocity[d]);
          }
          weight.push_back(p.weight);
          gamma.push_back(p.gamma);
        }
      }
    }

    writer.Put(var_loc, position.data());
    writer.Put(var_vel, velocity.data());
    writer.Put(var_w, weight.data());
    writer.Put(var_gamma, gamma.data());

    writer.EndStep();
    writer.Close();
  }

  ParticleMetric::ParticleMetric(const group_t* g_,
                                 adios2::IO&& io_,
                                 const std::size_t ncx,
                                 const std::size_t ncy,
                                 const std::size_t ncz)
  : io(io_),
    group(g_),
    var_density(io.DefineVariable<compute_t>("Density", {ncx, ncy, ncz}, {0, 0, 0}, {ncx, ncy, ncz}, adios2::ConstantDims)),
    var_temp(io.DefineVariable<compute_t>("Temperature", {ncx, ncy, ncz}, {0, 0, 0}, {ncx, ncy, ncz}, adios2::ConstantDims)),
    density(g_->cells.size()),
    T_avg(g_->cells.size()),
    KE_total(g_->cells.size())
  {}

  void ParticleMetric::update_metrics() {
    using Chunk = particles::ParticleChunk;
    constexpr auto V_cell_inv = 1.0_fp / (dx * dy * dz);
    constexpr auto temp_coef = 2.0_fp / (3.0_fp * constants::q_e<compute_t>);

    const auto mc2 = group->mass * constants::c_sqr<compute_t>;

    std::ranges::fill(density, 0.0_fp);
    std::ranges::fill(T_avg, 0.0_fp);
    std::ranges::fill(KE_total, 0.0_fp);

    for (const auto& [chunks, idxs] : group->cells) {
      if (chunks.empty()) { continue; }
      const auto cid = get_cid(idxs[0], idxs[1], idxs[2]);
      for (const auto& [particles, active, moves] : chunks) {
        for (std::size_t pid = 0; pid < Chunk::n_particles; pid++) {
          if (!active.test(pid)) { continue; }
          const auto& p = particles[pid];
          density[cid] += p.weight;
          KE_total[cid] += p.weight * mc2 * (p.gamma - 1.0);
        }
      }
    }

    for (std::size_t i = 0; i < T_avg.size(); i++) {
      if (density[i] == 0.0_fp) { continue; }
      T_avg[i] = temp_coef * KE_total[i] / density[i];
    }

    for (auto& x: density) {
      x *= V_cell_inv;
    }
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

  void Metrics::write(const std::size_t step) const {
    static constexpr size_t padding = 10;
    std::string count_padded = std::to_string(step);
    count_padded.insert(count_padded.begin(), padding - count_padded.length(), '0');
    count_padded += file_ext; // default extension is .bp

    for (const auto& m: metrics) {
      m->write(data_dir, count_padded);
    }
  }
}

