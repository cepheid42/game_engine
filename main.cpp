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

  EM1D::advance();


  // EM2D::advance();

  // auto start = std::chrono::high_resolution_clock::now();
  // auto stop = std::chrono::high_resolution_clock::now() - start;
  // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop);
  // std::cout << "Execution time: " << duration.count() << " microseconds" << std::endl;
  return 0;
}
