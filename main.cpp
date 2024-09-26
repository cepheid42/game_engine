#include <fstream>
#include <cmath>
#include <iostream>
#include <chrono>

#define DEBUG

#include "electromagnetics/em_solver.h"

template<typename T>
void print_type() {
  std::cout << __PRETTY_FUNCTION__ << std::endl;
}

int main() {
  using emdata_t = em_data_tm<tf::types::Array2D<double>>;

  const size_t nx = 10;

  emdata_t em{nx, nx};

  EMTM::advance(em);


  // EM2D::advance();

  // auto start = std::chrono::high_resolution_clock::now();
  // auto stop = std::chrono::high_resolution_clock::now() - start;
  // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop);
  // std::cout << "Execution time: " << duration.count() << " microseconds" << std::endl;
  return 0;
}
