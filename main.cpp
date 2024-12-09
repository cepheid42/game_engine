#include <fstream>
#include <cmath>
#include <iostream>
#include <chrono>
#include <string>

// #define DBG_MACRO_DISABLE
#define NTHREADS 32
#define NTHREADS_BC 16

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
requires is_empty_field<Array, EmptyArray<typename Array::value_t, Array::dimension_t::value>>
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

fp_t ricker(const fp_t q) {
  constexpr auto Np = 20.0;
  constexpr auto Md = 2.5;

  const auto alpha = (M_PI * (cfl * q / Np - Md)) * (M_PI * (cfl * q / Np - Md));

  return (1.0 - 2.0 * alpha) * std::exp(-alpha);
  // return std::exp(-sqr((q - 30.0) / 10.0));
}

int main() {
  const auto start = std::chrono::high_resolution_clock::now();

  constexpr size_t nx = 100u + 2 * nPml + 2 * nHalo;
  constexpr size_t ny = 100u + 2 * nPml + 2 * nHalo;
  // constexpr size_t nz = 100u + 2 * nPml + 2 * nHalo;
  constexpr size_t nt = 1u;

  // emdata_t<double> em{nx, cfl};
  // bcdata_t<double> bc{em};

  emdata_t<double> em{nx, ny, cfl};
  bcdata_t<double> bc{em};

  // emdata_t<double> em{nx, ny, nz, cfl};
  // bcdata_t<double> bc{em};

  constexpr auto save_step = 4;
  size_t filecount = 0;
  for (size_t n = 0; n < nt; n++) {

    EMSolver<fp_t>::advance(em, bc);

    const auto rsrc = ricker(static_cast<fp_t>(n));

    // std::cout << rsrc << std::endl;

    // em.Ez[nx / 2 - 20] += rsrc;

    // em.Ex(nx / 2 - 20, ny / 2 - 20) += rsrc;
    // em.Ey(nx / 2 - 20, ny / 2 - 20) += 3.0 * rsrc;
    em.Ez(nx / 2 - 20, ny / 2 - 20) += rsrc;

    // em.Ex(nx / 2, ny / 2, nz / 2 - 20) += rsrc;
    // em.Ey(nx / 2, ny / 2, nz / 2 - 20) += rsrc;
    // em.Ez(nx / 2, ny / 2, nz / 2 - 20) += rsrc;

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
//         }
//       }
      filecount++;
    }
  }

  const auto stop = std::chrono::high_resolution_clock::now() - start;
  const auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop);
  std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;
  return 0;
}
