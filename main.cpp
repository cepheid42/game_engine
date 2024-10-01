#include <fstream>
#include <cmath>
#include <iostream>
#include <chrono>

#define DEBUG

#include "electromagnetics/em_solver.h"
#include "electromagnetics/em_data.h"
#include "aydenstuff/array.h"

using array_t = tf::types::Array2D<double>;
using empty_t = EmptyArray2D<double>;
using emdata_t = em_data_tm<array_t>;

using EzTL = TypeList<array_t, array_t, array_t, array_t, array_t, array_t, array_t>;
using HxTL = TypeList<array_t, empty_t, array_t, empty_t, array_t, array_t, empty_t>;
using HyTL = TypeList<array_t, array_t, empty_t, empty_t, array_t, array_t, empty_t>;


using EIX = FieldIntegratorNull;
using EIY = FieldIntegratorNull;
using EIZ = FieldIntegrator2D<EzTL, Derivative::DX, Derivative::DY, false>;

using HIX = FieldIntegrator2D<HxTL, Derivative::NoOp, Derivative::DY, true>;
using HIY = FieldIntegrator2D<HyTL, Derivative::DX, Derivative::NoOp, true>;
using HIZ = FieldIntegratorNull;

using EMSolver = Electromagnetics<EIX, EIY, EIZ, HIX, HIY, HIZ>;

int main() {

  emdata_t emdata{5, 5};
  EMSolver::advance(emdata);

  // auto start = std::chrono::high_resolution_clock::now();
  // auto stop = std::chrono::high_resolution_clock::now() - start;
  // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop);
  // std::cout << "Execution time: " << duration.count() << " microseconds" << std::endl;
  return 0;
}
