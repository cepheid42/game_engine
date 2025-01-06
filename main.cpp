#include <fstream>
#include <cmath>
#include <iostream>
#include <chrono>
#include <string>

// #define DBG_MACRO_DISABLE

#include "electromagnetics/electromagnetics.h"

template<typename Array>
void to_csv(const Array& arr, const size_t step, const std::string& name) {
  std::string count_padded = std::to_string(step);
  count_padded.insert(count_padded.begin(), 6 - count_padded.length(), '0');

  std::string filename{"../data/" + name + "_" + count_padded + ".csv"};
  std::ofstream file;
  file.open(filename.c_str());

  if constexpr (std::is_same_v<Array, std::array<double, nPml>>) {
    for (const auto& e: arr) {
      file << e << ", ";
    }
  } else if constexpr (Array::dimension_t::value == 1) {
    for (size_t i = 0; i < arr.nx(); i++) {
      file << arr[i];
      if (i < arr.nx() - 1) {
        file << ", ";
      }
    }
    file << std::endl;
  } else if constexpr (Array::dimension_t::value == 2) {
    for (size_t i = 0; i < arr.nx(); i++) {
      for (size_t j = 0; j < arr.ny(); j++) {
        file << arr(i, j);
        if (j < arr.ny() - 1) {
          file << ", ";
        }
      }
      file << std::endl;
    }
  } else {
    for (size_t i = 0; i < arr.nx(); i++) {
      for (size_t j = 0; j < arr.ny(); j++) {
        for (size_t k = 0; k < arr.nz(); k++) {
          file << arr(i, j, k);
          if (k < arr.nz() - 1) {
            file << ", ";
          }
        }
        file << '\n';
      }
    }
  }

  file.close();
}
template<typename Array>
requires tf::electromagnetics::traits::instance_of_type<tf::types::EmptyArray, Array>
void print_array(const Array&) {
  std::cout << "Empty Array." << std::endl;
}

template<typename Array>
void print_array(const Array& arr) {
  if constexpr (Array::dimension_t::value == 1) {
    for (size_t i = 0; i < arr.nx(); i++) {
      std::cout << arr[i] << ", ";
    }
    std::cout << std::endl;
  } else {
    for (size_t i = 0; i < arr.nx(); i++) {
      for (size_t j = 0; j < arr.ny(); j++) {
        std::cout << arr(i, j);
        if (j < arr.ny() - 1) {
          std::cout << ", ";
        }
      }
      std::cout << std::endl;
    }
  }
}


int main() {
  const auto start = std::chrono::high_resolution_clock::now();

  constexpr size_t nx = 100u + 2 * nPml + 2 * nHalo;
  constexpr size_t ny = 100u + 2 * nPml + 2 * nHalo;
  // constexpr size_t nz = 100u + 2 * nPml + 2 * nHalo;
  constexpr size_t nt = 400u;

  constexpr double dx = 1.0 / 99.0;

  constexpr auto c0 = 299792458.0;
  constexpr auto dt = cfl * dx / c0;

  emdata_t<double> em{nx, ny, dt, dx};
  bcdata_t<double> bc{em, dt, dx};

  // constexpr size_t x0 = nPml + 5;
  // constexpr size_t x1 = nx - nPml - 5;
  // using tf::electromagnetics::sources::ContinuousSource;
  // using tf::electromagnetics::sources::GaussianSource;
  // using tf::electromagnetics::sources::GaussianBeam;
  // constexpr auto omega = 2.0 * M_PI * c0 / (40.0 * dx);
  // em.beams.push_back(
  //   std::make_unique<GaussianBeam<Array2D<fp_t>>>(
  //     GaussianBeam<Array2D<fp_t>>{
  //       &em.Ez,
  //       -2.75e13, // amp
  //       0.17, // waist
  //       omega, // freq
  //       {0.5, 0.5}, // waist position
  //       x0, // x0
  //       x0 + 1, // x1
  //       x0, // y0
  //       x1, // y1,
  //       dx
  //     }
  //   )
  // );
  //
  // em.beams[0]->t_srcs.push_back(std::make_unique<ContinuousSource<fp_t>>(omega, 0.0, 1.0e30, 0.0));
  //
  // constexpr auto width = 7.58e-9;
  // em.beams[0]->t_srcs.push_back(std::make_unique<GaussianSource<fp_t>>(width, 2.0, 2.0 * width));

  using tf::electromagnetics::sources::CurrentSource;
  using tf::electromagnetics::sources::SpatialSource;
  using tf::electromagnetics::sources::RickerSource;
  using temporal_vec = std::vector<std::unique_ptr<tf::electromagnetics::sources::TemporalSource<fp_t>>>;

  constexpr auto freq = c0 / (40.0 * dx);

  auto make_tvec = [&]()
  {
    temporal_vec tvec{};
    tvec.push_back(std::make_unique<RickerSource<fp_t>>(freq));
    return tvec;
  };

  em.srcs.push_back(
    std::make_unique<CurrentSource<Array2D<fp_t>>>(
      &em.Ez,
      SpatialSource<fp_t>(
        make_tvec(),
        50.0, 60, 61, 60, 61
      )
    )
  );


  // emdata_t<double> em{nx, ny, nz, dt, dx};
  // bcdata_t<double> bc{em, dt, dx};

  constexpr auto save_step = 4;
  size_t filecount = 0;
  for (size_t n = 0; n < nt; n++) {

    EMSolver<fp_t>::advance(static_cast<fp_t>(n) * dt, em, bc);

    if (n % save_step == 0) {
      std::cout << "Step " << n << std::endl;

// #pragma omp parallel num_threads(3)
//       {
// #pragma omp single
//         {
// #pragma omp task
//           to_csv(em.Ex, filecount, "Ex");
// #pragma omp task
//           to_csv(em.Ey, filecount, "Ey");
// #pragma omp task
           to_csv(em.Ez, filecount, "Ez");
// #pragma omp task
//           to_csv(em.Hx, filecount, "Hx");
// #pragma omp task
//           to_csv(em.Hy, filecount, "Hy");
// #pragma omp task
//           to_csv(em.Hz, filecount, "Hz");
//          }
//        }

      filecount++;
    }
  }

  const auto stop = std::chrono::high_resolution_clock::now() - start;
  const auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop);
  std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;
  return 0;
}
