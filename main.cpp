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
      file << arr[i] << ", ";
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

fp_t ricker(fp_t q) {
  constexpr auto Np = 20.0;
  constexpr auto Md = 1.0;

  const auto alpha = (M_PI * (cfl * q / Np - Md)) * (M_PI * (cfl * q / Np - Md));

  return (1.0 - 2.0 * alpha) * std::exp(-alpha);
  // return std::exp(-sqr((q - 30.0) / 10.0));
}

int main() {
  const auto start = std::chrono::high_resolution_clock::now();

  constexpr size_t nx = 100u + 2 * nPml + 2 * nHalo;
  constexpr size_t ny = 100u + 2 * nPml + 2 * nHalo;
  constexpr size_t nz = 100u + 2 * nPml + 2 * nHalo;
  constexpr size_t nt = 1000u;

  // emdata_t<double> em{nx, cfl};
  // bcdata_t<double> bc{em};

  // emdata_t<double> em{nx, ny, cfl};
  // bcdata_t<double> bc{em};

  emdata_t<double> em{nx, ny, nz, cfl};
  bcdata_t<double> bc{em};

  // constexpr auto eps0 = 8.854187812813e-12;
  // constexpr auto mu0 = 1.0 / 1.2566370621219e-6;

  constexpr auto save_step = 10;
  size_t filecount = 0;
  for (size_t n = 0; n < nt; n++) {

    EMSolver<fp_t>::advance(em, bc);

    // em.Ez[nx / 2] += ricker(static_cast<fp_t>(n));
    // em.Ex(nx / 2, ny / 2) += ricker(static_cast<fp_t>(n));
    // em.Ez(nx / 2, ny / 2) += ricker(static_cast<fp_t>(n));

    const auto rsrc = ricker(static_cast<fp_t>(n));
    // em.Ex(nx / 2, ny / 2, nz / 2) = rsrc;
    // em.Ey(nx / 2, ny / 2, nz / 2) = rsrc;
    em.Ez(nx / 2, ny / 2, nz / 2) = 1.0E4 * rsrc;

    // em.Hx(nx / 2, ny / 2, nz / 2) = rsrc / 1000.0;
    // em.Hy(nx / 2, ny / 2, nz / 2) = rsrc / 1000.0;
    // em.Hz(nx / 2, ny / 2, nz / 2) = rsrc / 1000.0;

    if (n % save_step == 0) {
      std::cout << "Step " << n << std::endl;
      to_csv(em.Ex, filecount, "Ex");
      to_csv(em.Ey, filecount, "Ey");
      to_csv(em.Ez, filecount, "Ez");
      to_csv(em.Hx, filecount, "Hx");
      to_csv(em.Hy, filecount, "Hy");
      to_csv(em.Hz, filecount, "Hz");
      filecount++;
    }
  }

  const auto stop = std::chrono::high_resolution_clock::now() - start;
  const auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop);
  std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;
  return 0;
}
