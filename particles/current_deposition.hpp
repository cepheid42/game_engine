#ifndef CURRENT_DEPOSITION_HPP
#define CURRENT_DEPOSITION_HPP

#include "math_utils.hpp"
#include "particles.hpp"
#include "em_data.hpp"
#include "bc_data.hpp"
#include "morton.hpp"

#include <array>
#include <cmath>
#include <cassert>

namespace tf::particles {
   struct TSC {
      static auto inner_shape(const compute_t x) {
         return 0.75_fp - math::SQR(std::abs(x));
      }

      static auto outer_shape(const compute_t x) {
         return 0.5_fp * math::SQR(1.5_fp - std::abs(x));
      }

      static auto operator()(const compute_t x) {
         if (std::abs(x) < 0.5_fp) {
            return inner_shape(x);
         }
         return outer_shape(x);
      }
   };

   template<bool isJ>
   struct ShapeFunctorJIT {
      compute_t particle_start{};
      compute_t particle_end{};
      TSC shape{};

      explicit ShapeFunctorJIT(const compute_t start_, const compute_t end_)
      : particle_start(start_), particle_end(end_)
      {}

      [[nodiscard]] auto S0(const int gridPoint) const {
         if constexpr (isJ) {
            return 1.0_fp;
         } else {
            return shape(gridPoint - particle_start);
         }
      }

      [[nodiscard]] auto S1(const int gridPoint) const {
         if constexpr (isJ) {
            return 1.0_fp;
         } else {
            return shape(gridPoint - particle_end);
         }
      }

      [[nodiscard]] auto DS(const int gridPoint) const {
         return S1(gridPoint) - S0(gridPoint);
      }
   };

   template<bool isJ>
   struct ShapeFunctorCached {
      std::array<compute_t, 3> shape0{};
      std::array<compute_t, 3> shape1{};
      TSC shape{};

      explicit ShapeFunctorCached(const compute_t x0, const compute_t x1)
      : shape0{TSC::outer_shape(x0), TSC::inner_shape(x0), 1.0_fp - (TSC::outer_shape(x0) + TSC::inner_shape(x0))},
        shape1{TSC::outer_shape(x1), TSC::inner_shape(x1), 1.0_fp - (TSC::outer_shape(x1) + TSC::inner_shape(x1))}
      {}

      [[nodiscard]] auto S0(const int gridPoint) const {
         if constexpr (isJ) {
            return 1.0;
         } else {
            return shape0[gridPoint + 1];
         }
      }

      [[nodiscard]] auto S1(const int gridPoint) const {
         if constexpr (isJ) {
            return 1.0;
         } else {
            return shape1[gridPoint + 1];
         }
      }

      [[nodiscard]] auto DS(const int gridPoint) const {
         return S1(gridPoint) - S0(gridPoint);
      }
   };

   template<int newX, int newY, int newZ>
   constexpr auto rotateOrigin(const auto& vec) {
      return decltype(vec){vec[newX], vec[newY], vec[newZ]};
   }

   template<int D>
   void applyJ(auto& J, const auto& acc, const std::array<std::size_t, 3>& cids, std::array<int, 3>&& ids) {
      const auto& [ci, cj, ck] = cids;

      if constexpr (D == 0) {
         auto [i, j, k] = rotateOrigin<2, 0, 1>(ids); // unrotates (1,2,0)->(0,1,2)
         if (j == -1) {
            j = 0;
         }
         // #pragma omp atomic update
         J(i + ci, j + cj, k + ck) += acc;
      }
      if constexpr (D == 1) {
         auto [i, j, k] = rotateOrigin<1, 2, 0>(ids); // unrotates (2,0,1)->(0,1,2)
         if (j == -1) {
            j = 0;
         }
         // #pragma omp atomic update
         J(i + ci, j + cj, k + ck) += acc;
      }
      if constexpr (D == 2) {
         auto [i, j, k] = ids;
         if (j == -1) {
            j = 0;
         }
         // #pragma omp atomic update
         J(i + ci, j + cj, k + ck) += acc;
      }
   }

   struct CurrentDeposition {
      using emdata_t = electromagnetics::EMData;
      using group_t = ParticleGroup;
      using array_t = std::array<std::size_t, 3>;
      using EMFace = electromagnetics::EMFace;
      using EMSide = electromagnetics::EMSide;

      template<int D>
      static void updateJ(auto& J, const auto qA,
                          auto& p0, auto& p1,
                          const std::array<std::size_t, 3>& idxs)
      {
         if (p0[D] == p1[D]) { return; }

         int i0, i1, j0, j1, k0, k1;
         if constexpr (D == 0) {
            p0 = rotateOrigin<1, 2, 0>(p0);
            p1 = rotateOrigin<1, 2, 0>(p1);
            i0 = 0;
            i1 = 2;
            j0 = -1;
            j1 = 2;
            k0 = -1;
            k1 = 1;
         } else {
            i0 = -1;
            i1 = 2;
            j0 = 0;
            j1 = 2;
            k0 = -1;
            k1 = 1;
         }

         const ShapeFunctorJIT<D==0>    shapeI(p0[0], p1[0]);
         const ShapeFunctorCached<D==2> shapeJ(p0[1], p1[1]); // template makes j-comp return 1.0
         const ShapeFunctorCached<false> shapeK(p0[2], p1[2]);

         for (int i = i0; i < i1; i++) {
            const auto s0i = shapeI.S0(i);
            const auto dsi = shapeI.S1(i) - s0i;
            for (int j = j0; j < j1; j++) {
               const auto s0j = shapeJ.S0(j);
               const auto dsj = shapeJ.S1(j) - s0j;
               const auto temp = qA * (s0i * s0j + 0.5_fp * (dsi * s0j + s0i * dsj) + (1.0_fp / 3.0_fp) * dsj * dsi);
               auto acc = 0.0_fp;
               for (int k = k0; k < k1; k++) {
                  acc += shapeK.DS(k) * temp;
                  applyJ<D>(J, acc, idxs, {i, j, k});
               }
            }
         }
      } // end updateJ()

      static void updateJy(auto& Jy, const auto qA,
                           auto& p0, auto& p1,
                           const std::array<std::size_t, 3>& idxs)
      {
         if (qA == 0.0) { return; }

         // p0 = rotateOrigin<2, 0, 1>(p0);
         // p1 = rotateOrigin<2, 0, 1>(p1);

         const auto shapeI = ShapeFunctorCached<false>(p0[0], p1[0]);
         const auto shapeK = ShapeFunctorCached<false>(p0[2], p1[2]);

         for (int i = -1; i < 2; i++) {
            const auto s0i = shapeI.S0(i);
            const auto dsi = shapeI.S1(i) - s0i;
            for (int k = -1; k < 2; k++) {
               const auto s0k = shapeK.S0(k);
               const auto dsk = shapeK.S1(k) - s0k;

               const auto W = s0i * s0k + 0.5_fp * (dsi * s0k + s0i * dsk) + (1.0_fp / 3.0_fp) * dsi * dsk;
               applyJ<1>(Jy, W * qA, idxs, {i, k, 0});
            }
         }
      }


      static void update_particle(const Particle& p, emdata_t& emdata, const compute_t charge) {
         constexpr auto dtAyz = 1.0_fp / (Ayz * dt);
         // constexpr auto dtAxz = 1.0_fp / (Axz * dt);
         constexpr auto dtAxy = 1.0_fp / (Axy * dt);

         auto relayPoint = [](const auto& a, const auto& b, const auto& x) { return (a == b) ? x : (a + b) / 2.0_fp; };
         const auto idxs = morton_decode(p.code);

         vec3<compute_t> i1{}, i2{}, p0{}, p1{};
         for (std::size_t d = 0; d < 3; ++d) {
            i1[d] = std::floor(p.old_location[d] + 0.5_fp);
            i2[d] = std::floor(p.location[d] + 0.5_fp);
            p0[d] = p.old_location[d] - i1[d];
            p1[d] = relayPoint(i1[d], i1[d], p.location[d]) - i1[d];
         }

         updateJ<0>(emdata.Jx, p.weight * charge * dtAyz, p0, p1, idxs);
         updateJy(emdata.Jy, p.velocity[1] * p.weight * charge, p0, p1, idxs); // 2D: y is done different
         updateJ<2>(emdata.Jz, p.weight * charge * dtAxy, p0, p1, idxs);

         if (i1 != i2) {
            vec3 p2 = p1 - i2;
            vec3 p3 = p.location - i2;

            updateJ<0>(emdata.Jx, p.weight * charge * dtAyz, p2, p3, idxs);
            updateJ<2>(emdata.Jz, p.weight * charge * dtAxy, p2, p3, idxs);
         }
      }

      static void operator()(const group_t& g, emdata_t& emdata) {
         // #pragma omp parallel for num_threads(nThreads)
         for (std::size_t pid = 0; pid < g.num_particles(); pid++) {
            update_particle(g.particles[pid], emdata, g.charge);
         }
      }
   }; // end struct CurrentDeposition
} // end namepsace tf::particles
#endif //CURRENT_DEPOSITION_HPP
