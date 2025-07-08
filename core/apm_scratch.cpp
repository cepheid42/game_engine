
// Universal includes
#include <iostream>
#include <concepts>

// Specialized includes
#include <algorithm>

#include "../libs/params/apm.select"
#include "../libs/common/include/array.h"
#include "../libs/common/include/constants.h"
#include "../libs/common/include/utilities/control.h"
#include "../libs/pic/include/pic_histogram.h"
#include "../libs/pic/include/pic_interpolation.h"
#include "../libs/pic/include/pic_particle.h"
#include "../libs/common/include/utilities/utilities.h"

// ===== Configuration =====
// =========================
namespace tf::options::apm
{
  constexpr auto  maximum_iterations = tf::select::apm::APM_MAX_ITER_SELECT;
  constexpr auto       minimum_count = tf::select::apm::APM_MIN_COUNT_SELECT;
  constexpr auto relative_low_weight = static_cast<double>(tf::select::apm::APM_REL_LOW_WEIGHT_SELECT);
  constexpr auto       buffer_factor = static_cast<double>(tf::select::apm::APM_BUFFER_FACTOR_SELECT);
  constexpr auto    target_threshold = static_cast<double>(tf::select::apm::APM_TARGET_THRESHOLD_SELECT);
  constexpr auto              margin = static_cast<double>(tf::select::apm::APM_MARGIN_SELECT);

  constexpr auto       subcells_per_cell = tf::select::apm::APM_SUBCELLS_PER_CELL_SELECT;
  constexpr auto total_subcells_per_cell = subcells_per_cell[0] * subcells_per_cell[1];
  //
} // end namespace tf::options::apm


// ----- APM Vector -----
template <typename T>
using APMVector = tf::types::Array1D<T>;

template <typename T>
auto& operator<<(std::ostream& os, APMVector<T> vec)
{
  auto nx = vec.nx();

  for (std::size_t i = 0; i < nx; i++) {
    os << vec(i);
    if (i != nx-1) os << " ";
  }
  os << std::endl;

  return os;
}

// ----- APM Matrix -----
template <typename T>
using APMMatrix = tf::types::Array2D<T>;

template <typename T>
auto& operator<<(std::ostream& os, APMMatrix<T> mat)
{
  auto nx = mat.nx();
  auto ny = mat.nz();

  for (std::size_t i = 0; i < nx; i++) {
    for (std::size_t j = 0; j < ny; j++) {
      os << mat(i, j);
      if (j != ny-1) os << " ";
    }
    os << std::endl;
  }

  return os;
}

// ----- APM Utils -----
template <typename T>
void swap_row(std::size_t a, std::size_t b, APMMatrix<T>& mat)
{
  auto aStart = mat.begin() + a * mat.nz();
  auto bStart = mat.begin() + b * mat.nz();
  std::swap_ranges(aStart, aStart + mat.nz(), bStart);
  //
} // end swap_row()

template <typename T>
void matmul(APMMatrix<T>& mat, APMVector<T>& vec, APMVector<T>& result)
{
  assert(vec.nx() == mat.nz());
  result = APMVector<T>{mat.nx()};

  auto m = mat.nx();
  auto n = mat.nz();
  for (size_t i = 0; i < m; i++) {
    result(i) = static_cast<T>(0.0);
    for (size_t j = 0; j < n; j++) {
      result(i) += mat(i, j) * vec(j);
    }
  }
  //
} // end matmul()

template <typename T>
bool lu_decomp(APMMatrix<T>& A, APMVector<std::size_t>& perm)
{
  static constexpr bool ERROR = true;
  static constexpr bool SUCCESS = false;

  // Initialize permutation vector
  auto s = perm.nx();
  for (std::size_t i = 0; i < s; i++) perm(i) = i;

  for (std::size_t j = 0; j < s-1; j++) {

    // Find index (k >= j) of max pivot
    auto pivot = j;
    for (std::size_t k = j; k < s; k++) {
      bool new_pivot = std::abs(A(k, j)) > std::abs(A(pivot, j));
      pivot = new_pivot ? k : pivot;
    } // end for(k)

    if (pivot != j) {
      swap_row(j, pivot, A);
      std::swap(perm(j), perm(pivot));
    }

    if (std::abs(A(j, j)) < 1.0e-6) { return ERROR; } // return error if matrix is nearly singular

    for (std::size_t k = j+1; k < s; k++) {
      A(k, j) = A(k, j) / A(j, j); // Jth column of L
      for (std::size_t m = j+1; m < s; m++) {
        A(k, m) = A(k, m) - ( A(k, j) * A(j, m) ); // Kth row of U
      } // end for(m)
    } // end for(k)

  } // end for(j)

  // Final singularity test
  for (std::size_t i = 0; i < s; i++) {
    if (std::abs(A(i, i)) < 1.0e-6) return ERROR;
  } // end for(i)

  return SUCCESS;
  //
} // end lu_decomp()

template <typename T>
void lu_solve(APMMatrix<T>& A, APMVector<std::size_t>& perm, APMVector<T>& b)
{
  auto s = p.size();

  auto temp = APMVector<T>(s);

  for (std::size_t i = 0; i < s; i++) { temp[i] = b[p[i]]; }
  b = std::move(temp);
  
  for (std::size_t i = 0; i < s; i++) {
    auto sum = 0.0;
    for (std::size_t j = 0; j < i; j++) {
      sum += A(i, j) * b[j];
    } // end for(j)
    b[i] -= sum;
  } // end for(i)

  b[s-1] = b[s-1]/A(s-1,s-1);
  for (std::size_t i = s-2; i >= 0; i--) {
    auto sum = 0.0;
    for (std::size_t j = s-1; j > i; j--) {
      sum += A(i, j) * b[j];
    } // end for(j)
    b[i] = (b[i] - sum) / A(i, i);
  } // end for(i)
  //
} // end lu_solve()

bool nr_solve(
  auto func,
  auto deriv,
  double& x,
  double min_value = -std::numeric_limits<double>::max(),
  double max_value = std::numeric_limits<double>::max(),
  double tolerance = 1e-6,
  const size_t max_inter = 10,
  double relax = 1.0
) {
  double x1 = 0.0;

  for (std::size_t iter = 0; iter < max_iter; iter++) {
    x1 = (1.0 - relax) * x + relax * (x - func(x) / deriv(x));

    const bool out_of_bounds = (x1 < min_value) || (x1 > max_value);
    double rel_change = std::abs((x1 - x) / x);

    x = x1;

    if (out_of_bounds) break;

    if (rel_change < tolerance) return true;
  }

  return false;
  //
} // end nr_solve()

// <AJK> all three of the following functions should be updated to use pre-computed nodal shapes
template <typename T, tf::concepts::Particle P, tf::concepts::Mesh M>
auto generate_apm_matrix(std::vector<P>& particles, M& mesh)
{
  using MyShape = tf::types::shapes::CIC<T>;
  APMMatrix<T> A{8, 8};

  for (auto& particle : particles)
  {
    auto index = tf::utilities::mesh::findIndex(particle.location, mesh);
    auto loc_norm = tf::utilities::interpolation::calculateNormalizedLocation(particle.location, mesh, index);

    auto x_shape = MyShape::shapeArray(loc_norm[0]);
    auto y_shape = MyShape::shapeArray(loc_norm[1]);
    auto z_shape = MyShape::shapeArray(loc_norm[2]);

    std::array<T, 8> node_shapes = {
      x_shape[0] * y_shape[0] * z_shape[0],
      x_shape[0] * y_shape[0] * z_shape[1],
      x_shape[0] * y_shape[1] * z_shape[0],
      x_shape[0] * y_shape[1] * z_shape[1],
      x_shape[1] * y_shape[0] * z_shape[0],
      x_shape[1] * y_shape[0] * z_shape[1],
      x_shape[1] * y_shape[1] * z_shape[0],
      x_shape[1] * y_shape[1] * z_shape[1]
    };

    for (std::size_t i = 0; i < 8; i++) {
      for (std::size_t j = 0; j < 8; j++) {
        A(i, j) += node_shapes[i] * node_shapes[j];
      } // end for(i)
    } // end for(j)
  } // end for(particle)

  return A;
  //
} // end generate_apm_matrix()

template <typename T, tf::concepts::Particle P, tf::concepts::Mesh M>
auto deposit_quantity(std::vector<P>& particles, M& mesh, auto& getter)
{
  std::array<T, 8> quantity;
  quantity.fill(0.0);

  for (auto& particle : particles) {
    auto index = tf::utilities::mesh::findIndex(particle.location, mesh);
    auto loc_norm = tf::utilities::interpolation::calculateNormalizedLocation(particle.location, mesh, index);

    auto x_shape = MyShape::shapeArray(loc_norm[0]);
    auto y_shape = MyShape::shapeArray(loc_norm[1]);
    auto z_shape = MyShape::shapeArray(loc_norm[2]);

    std::array<T, 8> node_shapes = {
      x_shape[0] * y_shape[0] * z_shape[0],
      x_shape[0] * y_shape[0] * z_shape[1],
      x_shape[0] * y_shape[1] * z_shape[0],
      x_shape[0] * y_shape[1] * z_shape[1],
      x_shape[1] * y_shape[0] * z_shape[0],
      x_shape[1] * y_shape[0] * z_shape[1],
      x_shape[1] * y_shape[1] * z_shape[0],
      x_shape[1] * y_shape[1] * z_shape[1]
    };

    for (std::size_t i = 0; i < 8; i++) {
      quantity[i] += node_shapes[i] * getter(particle);
    } // end for(i)
  } // end for(particle)

  return quantity;
  //
};// end deposit_quantity()

template <typename T, tf::concepts::Particle P, tf::concepts::Mesh M>
auto deposit_alternate_quantity(std::vector<P>& particles, M& mesh, auto& getter)
{
  auto numParticles = particles.size();

  std::array<T, 8> quantity;
  quantity.fill(0.0);

  for (std::size_t pid = 0; pid < numParticles; pid++) {
    const auto& particle = particles[pid];

    auto index = tf::utilities::mesh::findIndex(particle.location, mesh);
    auto loc_norm = tf::utilities::interpolation::calculateNormalizedLocation(particle.location, mesh, index);

    auto x_shape = MyShape::shapeArray(loc_norm[0]);
    auto y_shape = MyShape::shapeArray(loc_norm[1]);
    auto z_shape = MyShape::shapeArray(loc_norm[2]);

    std::array<T, 8> node_shapes = {
      x_shape[0] * y_shape[0] * z_shape[0],
      x_shape[0] * y_shape[0] * z_shape[1],
      x_shape[0] * y_shape[1] * z_shape[0],
      x_shape[0] * y_shape[1] * z_shape[1],
      x_shape[1] * y_shape[0] * z_shape[0],
      x_shape[1] * y_shape[0] * z_shape[1],
      x_shape[1] * y_shape[1] * z_shape[0],
      x_shape[1] * y_shape[1] * z_shape[1]
    };

    for (std::size_t i = 0; i < 8; i++) {
      quantity[i] += node_shapes[i] * getter(pid);
    } // end for(i)
  } // end for(particle)

  return quantity;
  //
};// end deposit_quantity()

template <typename T, tf::concepts::Particle P, tf::concepts::Mesh M>
void gather_quantity(std::vector<P>& particles, M& mesh, std::array<T, 8> quantity, auto& setter)
{
  for (auto& particle : particles) {
    auto index = tf::utilities::mesh::findIndex(particle.location, mesh);
    auto loc_norm = tf::utilities::interpolation::calculateNormalizedLocation(particle.location, mesh, index);

    auto x_shape = MyShape::shapeArray(loc_norm[0]);
    auto y_shape = MyShape::shapeArray(loc_norm[1]);
    auto z_shape = MyShape::shapeArray(loc_norm[2]);

    std::array<T, 8> node_shapes = {
      x_shape[0] * y_shape[0] * z_shape[0],
      x_shape[0] * y_shape[0] * z_shape[1],
      x_shape[0] * y_shape[1] * z_shape[0],
      x_shape[0] * y_shape[1] * z_shape[1],
      x_shape[1] * y_shape[0] * z_shape[0],
      x_shape[1] * y_shape[0] * z_shape[1],
      x_shape[1] * y_shape[1] * z_shape[0],
      x_shape[1] * y_shape[1] * z_shape[1]
    };

    auto sum = T(0.0);
    for (std::size_t i = 0; i < 8; i++) {
      sum += quantity[i] * node_shapes[i];
    } // end for(i)

    setter(particle, sum);
  } // end for(particle)

  //
} // end gather_quantity()

template <typename T, tf::concepts::Particle P>
void create_distribution(std::vector<P>& particles, tf::utilities::apm::Histogram<T>& histogram, auto& getter)
{
  // Weight & Mean Quantity
  T totalWeight = 0.0;
  T mean = 0.0;
  for (auto& particle : particles) {
    totalWeight += particle.weight;
    mean += particle.weight * getter(particle);
  }
  mean /= totalWeight;

  // Variance
  T variance = 0.0;
  for (auto& particle : particles) {
    variance += particle.weight * tf::utilities::math::SQR(getter(particle) - mean);
  }
  variance /= totalWeight;
  T deviation = std::sqrt(variance);

  histogram.set_min(mean - 3.0 * deviation);
  histogram.set_max(mean + 3.0 * deviation);

  for (auto& particle : particles) {
    histogram.accumulate(getter(particle), particle.weight);
  }

  histogram.update_cdf();
  //
} // end create_distribution()

// ===== APM Class =====
// =====================
enum APMError : int {SUCCESS = 0, WEIGHT_FAILURE = 1, ENERGY_FAILURE = 2};

template <tf::concepts::ParticleGroup G, tf::concepts::Mesh M>
class APM
{
  using particle_t = typename G::particle_t;

  explicit APM()
  : lowerThreshold((1.0 - tf::options::apm::target_threshold) * target_particles_per_cell),
    upperThreshold((1.0 + tf::options::apm::target_threshold) * target_particles_per_cell),
    rng{0.0, 1.0}
  {}

  ~APM() = default;

  void operator() (G& group, M& mesh, const time_t sim_time, const step_t step)
  {
    using tf::options::apm;

    if (interval.update_this_step(0.0, sim_time, step)) {
      size_t numSuccess = 0, numWeightFailure = 0, numEnergyFailure = 0; 

      // <AJK> these should probably be cached
      const auto lowerThreshold = target_particles_per_cell * (1.0 - target_threshold);
      const auto upperThreshold = target_particles_per_cell * (1.0 + target_threshold);
      
      for (std::size_t i = 0; i < mesh.nx(); i++) {
        for (std::size_t j = 0; j < mesh.ny(); j++) {
          for (std::size_t k = 0; k < mesh.nz(); k++) {
            std::array<std::size_t, 4> index{i, j, k, k + mesh.nz() * (j + (mesh.ny() * i))};
            const auto& ptrVec = group.particle_ptrs_per_cell[scid];
            const auto currentNumParticles = ptrVec.size();

            bool triggerAPM = (currentNumParticles < lowerThreshold or currentNumParticles > upperThreshold)
                            and currentNumParticles >= minimum_count;
            if (!triggerAPM) continue;

            std::vector<particle_t> old_particles{}, new_particles{};
            for (auto* ptr : ptrVec) old_particles.emplace_back(*ptr);

            auto apmErrorFlag = manage_cell(index, mesh, old_particles, new_particles);
            switch (apmErrorFlag) {
              case APMError::WEIGHT_FAILURE:
                numWeightFailure++;
                cull_bad_particles(new_particles);
                break;

              case APMError::ENERGY_FAILURE:
                numEnergyFailure++;
                cull_bad_particles(new_particles);
                break;

              case APMError::SUCCESS:
              default:
                numSuccess++;
                deactivate_particles(ptrVec);
                group.add_particles(new_particles);
                break;
            } // end switch(apmErrorFlag)
          } // end for(k)
        } // end for(j)
      } // end for(i)

      print_summary(numWeightFailure, numEnergyFailure, numSuccess);
    } // end if(update_this_step)
    //
  } // end operator()

  void deactivate_particles(std::vector<particle_t*> ptrs)
  {
    for (auto* ptr : ptrs) { ptr->active = false; }
    //
  } // end deactivate_particles()

  APMError manage_cell(std::array<std::size_t, 4> index, M& mesh, std::vector<particle_t>& old_particles, std::vector<particle_t>& new_particles, double mass)
  {
    APMError err = APMError::SUCCESS;

    // Step 0: Create new positions
    create_subcell_distributed_positions(index, mesh, old_particles, new_particles);

    // Step 1: Calculate new weights
    err = calculate_new_weights(index, mesh, old_particles, new_particles);

    size_t numWeightIterations = 1;
    while (err && (numWeightIterations < tf::options::apm::maximum_iterations)) {
      create_subcell_distributed_positions(index, mesh, old_particles, new_particles);

      err = calculate_new_weights(index, mesh, old_particles, new_particles);
      numWeightIterations++;
    }

    if (!err) {
      // Step 2: Calculate new momenta (average)
      calculate_new_avg_momentum(index, mesh, old_particles, new_particles, mass);

      // Step 3: Add thermal correction
      err = refine_energy_distribution(index, mesh, old_particles, new_particles, mass);
    }

    return err;
    //
  } // end manage_cell()

  void create_subcell_distributed_positions(
    const std::array<std::size_t, 4> index,
    const M& mesh,
    const std::vector<particle_t> old_particles,
    std::vector<particle_t> new_particles
  ) {
    const auto xc = mesh.x[index[0]];
    const auto dx = mesh.dx[index[0]];

    const auto yc = mesh.y[index[1]];
    const auto dy = mesh.dy[index[1]];

    const auto zc = mesh.z[index[2]];
    const auto dz = mesh.dz[index[2]];

    const auto oldNumParticles = old_particles.size();

    const float oldTotalWeight = 0.0
    for (const auto& particle : old_particles) { oldTotalWeight += particle.weight; }

    float targetWeight = oldTotalWeight / static_cast<float>(total_particles_per_cell);

    // Calculate subcell sizes
    auto sub_dx = dx / static_cast<float>(tf::options::apm::subcells_per_cell[0]);
    auto sub_dy = dy / static_cast<float>(tf::options::apm::subcells_per_cell[1]);
    auto sub_dz = dz / static_cast<float>(tf::options::apm::subcells_per_cell[2]);

    std::array<float, tf::options::apm::total_subcells_per_cell> totalSubCellWeight{0.0};
    for (const auto& particle : old_particles) {
      std::size_t si = std::min(static_cast<std::size_t>(std::floor((particle.location[0] - xc) / sub_dx)),
                                tf::options::apm::subcells_per_cell[0]-1);
      std::size_t sj = std::min(static_cast<std::size_t>(std::floor((particle.location[1] - yc) / sub_dy)),
                                tf::options::apm::subcells_per_cell[1]-1);
      std::size_t sk = std::min(static_cast<std::size_t>(std::floor((particle.location[2] - zc) / sub_dz)),
                                tf::options::apm::subcells_per_cell[2]-1);
      std::size_t sub_index = sk + tf::options::apm::subcells_per_cell[2] * (sj + (tf::options::apm::subcells_per_cell[1] * si));

      totalSubCellWeight[sub_index] += particle.weight;
    }

    std::array<std::size_t, tf::options::apm::total_subcells_per_cell> newSubCellCounts{0};
    for (std::size_t i = 0; i < tf::options::apm::total_subcells_per_cell; i++) {
      const std::size_t newCount = static_cast<size_t>(std::ciel(totalSubCellWeight[i] / targetWeight));
      newSubCellCounts[i] += newCount;
    }

    // Sample new positions
    for (std::size_t i = 0; i < tf::options::apm::subcells_per_cell[0]; i++) {
      for (std::size_t j = 0; j < tf::options::apm::subcells_per_cell[1]; j++) {
        for (std::size_t k = 0; k < tf::options::apm::subcells_per_cell[2]; k++) {
          auto xSubcell = xc + (static_cast<float>(i) * sub_dx);
          auto ySubcell = yc + (static_cast<float>(j) * sub_dy);
          auto zSubcell = zc + (static_cast<float>(k) * sub_dz);

          auto xNewParticle = xSubcell + (rng() * sub_dx);
          auto yNewParticle = ySubcell + (rng() * sub_dy);
          auto zNewParticle = zSubcell + (rng() * sub_dz);

          for (std::size_t n = 0; n < newSubCellCounts[0]; n++) {
            new_particles.emplace_back(
              location,
              old_location,
              momentum,
              gamma,
              weight
            )
          } // end for(n)
        } // end for(k)
      } // end for(j)
    } // end for(i)

    //
  } // end create_subcell_distributed_positions()

  bool calculate_new_weights(
    const std::array<std::size_t, 4> index,
    const M& mesh,
    const std::vector<particle_t> old_particles,
    std::vector<particle_t> new_particles
  )
  {
    constexpr auto get_weight = [](particle_t& particle) { return particle.weight; };
    auto q = deposit_quantity<double>(old_particles, mesh, get_weight);
    auto A = generate_apm_matrix<double>(new_particles, mesh);

    APMVector<std::size_t> p;
    if (lu_decomp(A, p)) { return true; }

    auto set_weight = [](particle_t& particle, double value) { particle.weight = value; }
    gather_quantity<double>(new_particles, mesh, q, set_weight);

    auto oldTotalWeight = 0.0;
    for (auto& particle : old_particles) { oldTotalWeight += particle.weight; }
    auto oldAvgWeight = oldTotalWeight / old_particles.size();

    auto newTotalWeight = 0.0;
    for (auto& particle : new_particles) {
      auto badWeightFlag = particle.weight <= (kRelativeLowWeight * oldAvgWeight);
      if (badWeightFlag) { return true; }

      newTotalWeight += particle.weight;
    }

    return false;
    //
  } // end calculate_new_weights()

  void calculate_new_avg_momentum(
    const std::array<std::size_t, 4> index,
    const M& mesh,
    const std::vector<particle_t> old_particles,
    std::vector<particle_t> new_particles,
    double mass
  )
  {
    // I NEED MASS FROM SOMEWHERE

    auto get_weighted_px = [&](particle_t& particle) { return particle.weight * particle.momentum[0] / mass / constants::c; };
    auto get_weighted_py = [&](particle_t& particle) { return particle.weight * particle.momentum[1] / mass / constants::c; };
    auto get_weighted_pz = [&](particle_t& particle) { return particle.weight * particle.momentum[2] / mass / constants::c; };
    auto set_avg_px = [&](particle_t& particle, double value) { particle.momentum[0] = value * mass * constants::c / particle.weight; };
    auto set_avg_py = [&](particle_t& particle, double value) { particle.momentum[1] = value * mass * constants::c / particle.weight; };
    auto set_avg_pz = [&](particle_t& particle, double value) { particle.momentum[2] = value * mass * constants::c / particle.weight; };

    auto qx = deposit_quantity<double>(old_particles, mesh, get_weighted_px);
    auto qy = deposit_quantity<double>(old_particles, mesh, get_weighted_py);
    auto qz = deposit_quantity<double>(old_particles, mesh, get_weighted_pz);
    auto A = generate_apm_matrix<double>(new_particles, mesh);

    APMVector<std::size_t> p;
    lu_decomp(A, p);
    lu_solve(A, p, qx);
    lu_solve(A, p, qy);
    lu_solve(A, p, qz);

    gather_quantity(new_particles, mesh, qx, set_avg_px);
    gather_quantity(new_particles, mesh, qy, set_avg_py);
    gather_quantity(new_particles, mesh, qz, set_avg_pz);
  }

  bool refine_energy_distribution(
    const std::array<std::size_t, 4> index,
    const M& mesh,
    const std::vector<particle_t> old_particles,
    std::vector<particle_t> new_particles,
    double mass
  )
  {
    // I NEED MASS FROM SOMEWHERE

    auto nBins = std::floor((old_particles.size() / 5) + 1);

    auto get_px_norm = [&](particle_t& particle) { return particle.momentum[0] / tf::constants::c / mass; };
    auto get_py_norm = [&](particle_t& particle) { return particle.momentum[1] / tf::constants::c / mass; };
    auto get_pz_norm = [&](particle_t& particle) { return particle.momentum[2] / tf::constants::c / mass; };
    auto get_p_norm = [&](particle_t& particle) { return particle.momentum / tf::constants::c / mass; };
    auto apply_px_correction = [&](particle_t& particle, double value) { return particle.mometum[0] += value * mass * tf::constants::c; };
    auto apply_py_correction = [&](particle_t& particle, double value) { return particle.mometum[1] += value * mass * tf::constants::c; };
    auto apply_pz_correction = [&](particle_t& particle, double value) { return particle.mometum[2] += value * mass * tf::constants::c; };

    auto px_dist = tf::utilities::apm::Histogram<double>(nBins);
    create_distribution(old_particles, px_dist, get_px_norm);

    auto py_dist = tf::utilities::apm::Histogram<double>(nBins);
    create_distribution(old_particles, py_dist, get_py_norm);

    auto pz_dist = tf::utilities::apm::Histogram<double>(nBins);
    create_distribution(old_particles, pz_dist, get_pz_norm);

    std::vector<double> px_accum(new_particles.size());
    std::vector<double> py_accum(new_particles.size());
    std::vector<double> pz_accum(new_particles.size());

    std::vector<double> px_rand(new_particles.size());
    std::vector<double> py_rand(new_particles.size());
    std::vector<double> pz_rand(new_particles.size());

    std::vector<double> c1(new_particles.size());
    std::vector<double> c2(new_particles.size());
    std::vector<double> c3(new_particles.size());

    bool isConverged = false;
    for (std::size_t iter = 0; iter < kMaxIterations; iter++) {
      std::generate(px_rand.begin(), px_rand.end(), [&]() { return px_dist.sample_from_histogram(); });
      std::generate(py_rand.begin(), py_rand.end(), [&]() { return py_dist.sample_from_histogram(); });
      std::generate(pz_rand.begin(), pz_rand.end(), [&]() { return pz_dist.sample_from_histogram(); });

      for (std::size_t pid = 0; pid < new_particles.size(); pid++) {
        px_rand[pid] -= get_px_norm(new_particles[pid]);
        py_rand[pid] -= get_py_norm(new_particles[pid]);
        pz_rand[pid] -= get_pz_norm(new_particles[pid]);
      }

      auto get_px_rand = [&](std::size_t pid) { return -1.0 * px_rand[pid]; }
      auto qx = deposit_alternate_quantity(new_particles, mesh, get_px_rand);

      auto get_py_rand = [&](std::size_t pid) { return -1.0 * py_rand[pid]; }
      auto qy = deposit_alternate_quantity(new_particles, mesh, get_py_rand);

      auto get_pz_rand = [&](std::size_t pid) { return -1.0 * pz_rand[pid]; }
      auto qz = deposit_alternate_quantity(new_particles, mesh, get_pz_rand);

      auto A = generate_apm_matrix(new_particles, mesh);

      APMVector<std::size_t> p;
      lu_decomp(A, p);
      lu_solve(A, p, qx);
      lu_solve(A, p, qy);
      lu_solve(A, p, qz);

      gather_quantity(new_particles, mesh, qx, apply_px_correction);
      gather_quantity(new_particles, mesh, qy, apply_py_correction);
      gather_quantity(new_particles, mesh, qz, apply_pz_correction);

      auto lhs = 0.0;
      auto maxEnergy = 0.0;
      for (std::size_t pid = 0; pid < old_particles.size(); pid++) {
        auto& particle = old_particles[pid];
        auto pNorm = get_p_norm(particle);

        auto energy = std::sqrt(1.0 + pNorm.length_squared()) - 1.0;
        lhs += particle.weight * (energy + 1.0);
        maxEnergy = std::max(maxEnergy, energy);
      }
      maxEnergy *= 1.10; // TODO <AJK> this scaling factor needs to be set somewhere

      c1.clear();
      c2.clear();
      c3.clear();
      for (std::size_t pid = 0; pid < new_particles.size(); pid++) {
        auto& particle = new_particles[pid];
        auto pNorm = get_p_norm(particle);

        c1.push_back(1.0 + p_norm.length_squared());

        c2.push_back((pNorm[0] * px_rand[pid]
                    + pNorm[1] * py_rand[pid]
                    + pNorm[2] * pz_rand[pid]) / particle.weight);

        c3.push_back((tf::utilities::math::SQR(px_rand[pid] / particle.weight)
                    + tf::utilities::math::SQR(py_rand[pid] / particle.weight)
                    + tf::utilities::math::SQR(pz_rand[pid] / particle.weight)));
      }

      auto f = [&](double beta)
      {
        auto sum = 0.0;
        for (std::size_t pid = 0; pid < new_particles.size(); pid++) {
          auto& particle = new_particles[pid];

          auto arg = c1[pid]
                   + (2.0 * c2[pid] * beta)
                   + (c3[pid] * tf::utilities::math::SQR(beta));

          sum += particle.weightg * std::sqrt(arg);
        }

        return sum - lhs;
      }

      auto df = [&](double beta)
      {
        auto sum = 0.0;
        for (std::size_t pid = 0; pid < new_particles.size(); pid++) {
          auto& particle = new_particles[pid];

          auto arg = c1[pid]
                   + (2.0 * c2[pid] * beta)
                   + (c3[pid] * tf::utilities::math::SQR(beta));
          sum += particle.weight * (c2[pid] + c3[pid] * beta) / std::sqrt(arg);
        }

        return sum;
      }

      auto meanWeight = 0.0;
      for(std::size_t pid = 0.0; new_particles.size(); pid++) {
        meanWeight += new_particles[pid].weight;
      }
      meanWeight /= new_particles.size();

      auto beta = 1.0 * meanWeight;
      isConverged = nr_solve(f, df, beta, -5.0 * meanWeight, 5.0 * meanWeight, 1.0e-6, 10, 1.0);

      if (isConverged) {
        for (std::size_t pid = 0; pid < new_particles.size(); pid++) {
          auto& particle = new_particles[pid];
          auto pNorm = get_p_norm(particle);
          auto betaOverWeight = beta / particle.weight;

          px_accum[pid] = betaOverWeight * px_rand[pid];
          py_accum[pid] = betaOverWeight * py_rand[pid];
          pz_accum[pid] = betaOverWeight * pz_rand[pid];

          auto energy = std::sqrt(1.0
                      + tf::utilities::math::SQR(pNorm[0] + px_accum[pid])
                      + tf::utilities::math::SQR(pNorm[1] + py_accum[pid])
                      + tf::utilities::math::SQR(pNorm[2] + pz_accum[pid])) - 1.0;

          if (energy > maxEnergy) {
            isConverged = false;
            break;
          }
        } // end for(pid)
      } // end if(isConverged)

      if (isConverged) {
        for (std::size_t pid = 0; pid < new_particles.size(); pid++) {
          auto& particle = new_particles[pid];
          auto pNorm = get_p_norm(particle);

          pNorm[0] += px_accum[pid];
          pNorm[1] += py_accum[pid];
          pNorm[2] += pz_accum[pid];
          particle.gamma = std::sqrt(1.0 + pNorm.length_squared());

          particle.momentum[0] = pNorm[0] * mass * tf::constants::c;
          particle.momentum[1] = pNorm[1] * mass * tf::constants::c;
          particle.momentum[2] = pNorm[2] * mass * tf::constants::c;
        } // end for(pid)
        break;
      } // end if(isConverged)
    } // end for(iter)

    return !isConverged;
    //
  } // end refine_energy_distribution()

  // Data Members
  tf::utilities::math::uniform_rng rng;
  std::size_t lowerThreshold, upperThreshold;
  tf::types::Interval interval;
  //
};// end class APM

int main()
{

  
  return 0;
}