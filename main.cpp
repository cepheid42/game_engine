#include <fstream>
#include <cmath>
#include <iostream>
#include <chrono>

#define DEBUG

#include "electromagnetics/em_solver.h"



int main() {

  em_data_1d<double> emdata{};
  Electromagnetics::advance(emdata);

  // auto start = std::chrono::high_resolution_clock::now();
  // auto stop = std::chrono::high_resolution_clock::now() - start;
  // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop);
  // std::cout << "Execution time: " << duration.count() << " microseconds" << std::endl;
  return 0;
}
