#include <fstream>
#include <cmath>
#include <iostream>
#include <chrono>

#define DEBUG

// #include "electromagnetics/em_solver.h"
#include "electromagnetics/em_data.h"
// #include "aydenstuff/array.h"

// #include "electromagnetics/bread.h"

// #include "electromagnetics/old_emdata_v2.h"
// #include "electromagnetics/old_solver_v2.h"

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

// using EzTL = TypeList<array_t, array_t, array_t, array_t, array_t, array_t, array_t>;
// using HxTL = TypeList<array_t, empty_t, array_t, empty_t, array_t, array_t, empty_t>;
// using HyTL = TypeList<array_t, array_t, empty_t, empty_t, array_t, array_t, empty_t>;
//
//
// using EIX = FieldIntegratorNull;
// using EIY = FieldIntegratorNull;
// using EIZ = FieldIntegrator2D<EzTL, Derivative::DX, Derivative::DY, false>;
//
// using HIX = FieldIntegrator2D<HxTL, Derivative::NoOp, Derivative::DY, true>;
// using HIY = FieldIntegrator2D<HyTL, Derivative::DX, Derivative::NoOp, true>;
// using HIZ = FieldIntegratorNull;
//
// using EMSolver = Electromagnetics<EIX, EIY, EIZ, HIX, HIY, HIZ>;

int main() {

  // emdata_t<double> emdata1d{5};
  //
  // emdata1d.Ez[0] = 1.0;
  // DBG(emdata1d.Ez);

  DBG(sizeof(emdata_t<double>));
  DBG(sizeof(emdata_t<double>::ex_t));
  DBG(sizeof(emdata_t<double>::ey_t));
  DBG(sizeof(emdata_t<double>::));
  DBG(sizeof(emdata_t<double>));
  // emdata_t<double> emdata2d{5, 5};
  //
  // emdata2d.Ez(0, 0) = 1.0;
  // DBG(emdata2d.Ez);

  // print_type<emdata_t<double>>();
  //
  // print_type<emdata_t<double>::ex_t>();
  // print_type<emdata_t<double>::ey_t>();
  // print_type<emdata_t<double>::ez_t>();



  // auto start = std::chrono::high_resolution_clock::now();
  // auto stop = std::chrono::high_resolution_clock::now() - start;
  // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop);
  // std::cout << "Execution time: " << duration.count() << " microseconds" << std::endl;
  return 0;
}
