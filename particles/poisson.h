#ifndef POISSON_H
#define POISSON_H

#include "program_params.hpp"
#include "array.hpp"

namespace tf::particles {

struct PoissonSolver {
   // 1. Backstep velocity
   // 2. Update particle old_locations
   // 3. Calculate J for t=0
   // 4. Solve for A using ∇²A = -μ₀J
   // 5. Solve for B using B = ∇⨯A
   // 6. Calculate ρ
   // 7. Calculate ϕ using ∇²ϕ = -ρ/ε₀
   // 8. Calculate E using E = -∇ϕ

   static auto createPoissonMatrix() {
      Array3D<compute_t> P{};

      return P;
   }

   static void calculateVectorPotential() {}

   static void calculateB(auto& emdata, const auto& P) {
      Array3D<compute_t> Ax{}; // Vector Potential
      Array3D<compute_t> Ay{}; // Vector Potential
      Array3D<compute_t> Az{}; // Vector Potential

      // Step 4
      calculateVectorPotential(emdata.Jx, emdata.Jy, emdata.Jz, Ax, Ay, Az, P);
   }

   static auto calculateRho(const auto& p_groups) {
      Array3D<compute_t> rho{};

      // todo: There is a metric for calculating this. Maybe worth making it second order or something...

      return rho;
   }

   static auto calculatePhi() {
      Array3D<compute_t> phi{};

      // Step 6
      const auto rho = calculateRho();

      return phi;
   }

   static void calculateE(auto& emdata, const auto& P) {
      // Step 7
      const auto phi = calculatePhi();
   }

   static void operator()(auto& emdata, const auto& p_groups) {
      // todo: should p_groups be a vector? Do poisson one at a time?

      const auto P = createPoissonMatrix();

      // Assume steps 1-3 are done already

      // steps 4 + 5 can be merged, since A isn't used for anything else
      // Step 5
      calculateB(emdata);

      // steps 6+7+8 can be merged into one, since I don't need phi or rho for anything else

      // Step 8
      calculateE(emdata);
   }
}; // end struct PoissonSolver
} // end namespace tf::particles

#endif //POISSON_H
