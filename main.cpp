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
  count_padded.insert(count_padded.begin(), 4 - count_padded.length(), '0');

  std::string filename{"../data/" + name + "_" + count_padded + ".csv"};
  std::ofstream file;
  file.open(filename.c_str());

  if constexpr (Array::dimension_t::value == 1) {
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
      file << std::endl;
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
  constexpr size_t nx = 100u + 2 * dPML + 2 * nHalo;
  constexpr size_t ny = 100u + 2 * dPML + 2 * nHalo;
  // constexpr size_t nz = 10u + 2 * dPML + 2 * nHalo;
  // constexpr size_t nt = 400u;


  // emdata_t<double> em{nx, cfl};
  // bcdata_t<double> bc{em};

  emdata_t<double> em{nx, ny, cfl};
  bcdata_t<double> bc{em};

  EMSolver<fp_t>::advance(em, bc);

  // constexpr auto save_step = 4;
  // size_t filecount = 0;
  // for (size_t n = 0; n < nt; n++) {
  //   std::cout << "Step " << n << std::endl;
  //
  //   EMSolver<fp_t>::advance(em, bc);
  //
  //   em.Ez[nx / 2 + 20] += ricker(static_cast<fp_t>(n));
  //   // em.Ez(nx / 2 - 30, ny / 2 - 30) = ricker(static_cast<fp_t>(n));
  //   // em.Ez(nx / 2 - 30, ny / 2 - 30, nz / 2) = ricker(static_cast<fp_t>(n));
  //
  //   if (n % save_step == 0) {
  //     to_csv(em.Ez, filecount, "Ez");
  //     // to_csv(em.Hy, filecount, "Hy");
  //     filecount++;
  //   }
  // }

  // auto start = std::chrono::high_resolution_clock::now();
  // auto stop = std::chrono::high_resolution_clock::now() - start;
  // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop);
  // std::cout << "Execution time: " << duration.count() << " microseconds" << std::endl;
  return 0;
}
