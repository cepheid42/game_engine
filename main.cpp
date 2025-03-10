#include <iostream>
// #include <vector>

#include "program_params.hpp"
#include "em_solver.hpp"

int main() {
  std::cout << Nx << ", " << Ny << ", " << Nz << std::endl;
  tf::electromagnetics::EMSolver emsolver(Nx, Ny, Nz, cfl, dt);

  for (std::size_t n = 0; n < 400; n++) {
    emsolver.advance();
  }

  return 0;
}
