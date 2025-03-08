#include <fstream>
#include <cmath>
#include <iostream>
// #include <algorithm>
// #include <ranges>
// #include <sstream>
// #include <functional>
// #include <vector>

#define DEBUG

// #include "core/debug.h"
// #include "core/h5_wrapper.h"
// #include "core/logging.h"
// #include "core/maths.h"
// #include "core/memory.h"
// #include "core/profiling.h"
// #include "core/random.h"
// #include "core/timers.h"
// #include "core/typelist.h"
// #include "old_em/em_solver.h"

#include "expr_templates/expr_array.h"
#include "expr_templates/expr_scalar.h"
#include "expr_templates/expr_ops.h"

void print(auto& arr) {
  const auto shape = arr.shape();
  for (size_t i = 0; i < shape.x; i++) {
    for (size_t j = 0; j < shape.y; j++) {
      std::cout << arr(i, j);
      if (j < shape.y - 1) { std::cout << ", "; }
    }
    std::cout << '\n';
  }
  std::cout << std::endl;
}

void test_expr1d() {
  constexpr size_t nx = 50000;
  // constexpr size_t ny = 5;

  Array<double> x(nx), y(nx), z(nx);

  Array<int> idxs(3);
  idxs[0] = 0;
  idxs[1] = 2;
  idxs[2] = 4;

  for (std::size_t i = 0; i < nx; i++) {
    x[i] = static_cast<double>(i);
    y[i] = static_cast<double>(nx - i - 1);
  }


  for (std::size_t i = 0; i < nx; i++) {
    // std::cout << 2 << " * " << x[i] << " + " << y[i] << " = " << 2.0 * x[i] - y[i] << std::endl;
    // std::cout << x[i] << " + " << 2.0  << " = " << x[i] + 2.0 << std::endl;
    // std::cout << 2.0 << " - " << x[i] << " = " << 2.0 - x[i] << std::endl;
    // std::cout << x[i] << " - " << 2.0 << " = " << x[i] - 2.0 << std::endl;
  }

  // z[idxs] = 2.0 * x[idxs] - y[idxs];
  // z = 2.0 - x;
  // z = x - 2.0;

// #pragma omp parallel
  z = 2.0 * x - y;


  // for (std::size_t i = 0; i < nx; i++) {
  //   std::cout << z[i] << std::endl;
  // }
}

// void test_expr2d() {
//   constexpr std::size_t nx = 5;
//   constexpr std::size_t ny = 5;
//
//   Array<double> x(nx, ny), y(nx, ny), z(nx, nx);
//
//   Array<int> idxs(3);
//   idxs[0] = 0;
//   idxs[1] = 2;
//   idxs[2] = 4;
//
//   for (std::size_t i = 0; i < nx; i++) {
//     for (std::size_t j = 0; j < ny; j++) {
//       x(i, j) = static_cast<double>(i);
//       y(i, j) = static_cast<double>(j);
//     }
//   }
//
//   // for (std::size_t i = 0; i < nx; i++) {
//   //   for (std::size_t j = 0; j < ny; j++) {
//   //     std::cout << 2.0 << " * " << x(i, j) << " + " << y(i, j) << " = " << 2.0 * x(i, j) - y(i, j) << '\n';
//   //   }
//   // }
//
//   z = x - 2.0;
//
//   // print(x);
//   // print(y);
//   print(z);
// }

#include <chrono>
int main() {
  auto start = std::chrono::high_resolution_clock::now();
  test_expr1d();
  // test_expr2d();
  auto stop = std::chrono::high_resolution_clock::now() - start;

  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop);

  std::cout << "Execution time: " << duration.count() << " microseconds" << std::endl;
  return 0;
}
