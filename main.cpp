#include <fstream>
#include <cmath>
#include <iostream>
#include <chrono>

#define DEBUG

#include "electromagnetics/em_solver.h"
#include "electromagnetics/em_data.h"

using fp_t = double;

template<typename T>
void print_type() {
  std::cout << __PRETTY_FUNCTION__ << std::endl;
}

// using array_t = tf::types::Array2D<double>;
// using empty_t = EmptyArray2D<double>;
//
// // using efields =
//
// using emdata_t = test<ex_field<array_t>, ey_field<array_t>, ez_field<array_t>>;

// using emdata_t = em_data_tm<array_t>;

using emdata_fp = emdata_t<fp_t>;

using EzTL = TypeList<emdata_fp::ez_t, emdata_fp::hy_t, emdata_fp::hx_t, emdata_fp::ez_t, emdata_fp::ez_t, emdata_fp::ez_t, emdata_fp::ez_t>;
using HxTL = TypeList<emdata_fp::hx_t, emdata_fp::ey_t, emdata_fp::ez_t, emdata_fp::empty_t, emdata_fp::hx_t, emdata_fp::hx_t, emdata_fp::empty_t>;
using HyTL = TypeList<emdata_fp::hy_t, emdata_fp::ez_t, emdata_fp::ex_t, emdata_fp::empty_t, emdata_fp::hy_t, emdata_fp::hy_t, emdata_fp::empty_t>;


using EIX = FieldIntegratorNull<emdata_t<fp_t>>;
using EIY = FieldIntegratorNull<emdata_t<fp_t>>;
using EIZ = FieldIntegrator2D<EzTL, Derivative::DX, Derivative::DY, false>;

using HIX = FieldIntegrator2D<HxTL, Derivative::NoOp, Derivative::DY, true>;
using HIY = FieldIntegrator2D<HyTL, Derivative::DX, Derivative::NoOp, true>;
using HIZ = FieldIntegratorNull<emdata_t<fp_t>>;

using EMSolver = Electromagnetics<EIX, EIY, EIZ, HIX, HIY, HIZ>;

int main() {
  emdata_t<double> emdata2d{5, 5};

  emdata2d.Ez(0, 0) = 1.0;
  DBG(emdata2d.Ez);

  // auto start = std::chrono::high_resolution_clock::now();
  // auto stop = std::chrono::high_resolution_clock::now() - start;
  // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop);
  // std::cout << "Execution time: " << duration.count() << " microseconds" << std::endl;
  return 0;
}
