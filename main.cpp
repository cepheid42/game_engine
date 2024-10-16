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

int main() {
  // emdata_t<double> em{3};
  // using EMType = EM1D<emdata_t<double>>;

  emdata_t<double> em{3, 3};
  using EMType = EMTM<emdata_t<double>>;

  using EMSolver = Electromagnetics<
    TypeListAt<0, EMType>,
    TypeListAt<1, EMType>,
    TypeListAt<2, EMType>,
    TypeListAt<3, EMType>,
    TypeListAt<4, EMType>,
    TypeListAt<5, EMType>
  >;


  EMSolver::advance(em);


  // auto start = std::chrono::high_resolution_clock::now();
  // auto stop = std::chrono::high_resolution_clock::now() - start;
  // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop);
  // std::cout << "Execution time: " << duration.count() << " microseconds" << std::endl;
  return 0;
}
