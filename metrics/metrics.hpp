#ifndef METRICS_HPP
#define METRICS_HPP

#include "program_params.hpp"
#include "array.hpp"
#include "particles.hpp"

#include <print>
#include <unordered_map>
#include <memory>

#include <adios2.h>

namespace tf::metrics {
  namespace detail {
    struct MetricBase {
      virtual ~MetricBase() = default;
      virtual void write(const std::string&, const std::string&) = 0;
    };
  } // end namespace detail


  struct EMFieldsMetric final : detail::MetricBase {
    using pointer_t = Array3D<compute_t>*;
    using field_map = std::unordered_map<std::string, pointer_t>;

    struct FieldVariable {
      pointer_t field;
      adios2::Variable<compute_t> variable;
    };

    EMFieldsMetric(const field_map&, adios2::IO&&);
    void write(const std::string&, const std::string&) override;

    adios2::IO io;
    std::vector<FieldVariable> fields;
  };

  struct ParticleMetric final : detail::MetricBase {
    using group_t = particles::ParticleGroup;

    ParticleMetric(const group_t*, adios2::IO&&);
    void write(const std::string&, const std::string&) override;

    adios2::IO io;
    const group_t* group;
    adios2::Variable<compute_t> var_loc;
    adios2::Variable<double> var_vel;
    adios2::Variable<compute_t> var_w;
  };

  class Metrics {
    using metrics_vec = std::vector<std::unique_ptr<detail::MetricBase>>;

  public:
    explicit Metrics(std::string data_dir_) : data_dir(std::move(data_dir_)) {}
    void addMetric(std::unique_ptr<detail::MetricBase>&& m) { metrics.push_back(std::move(m)); }

    void write(std::size_t);

  public:
    adios2::ADIOS adios{};

  private:
    std::string file_ext{".bp"};
    std::string data_dir{};
    metrics_vec metrics{};
  };
}

#endif //METRICS_HPP
