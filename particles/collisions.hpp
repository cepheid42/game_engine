#ifndef COLLISIONS_HPP
#define COLLISIONS_HPP

#include "binary_channels.hpp"
#include "particles.hpp"
#include "constants.hpp"
#include "math_utils.hpp"
#include "rng_utils.h"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <ranges>
#include <span>
#include <vector>
// #include <print>

namespace tf::collisions
{

struct Collisions {
   using group_t = particles::ParticleGroup;
   using particle_vec = std::vector<particles::Particle>;
   group_t& g1; // electrons
   group_t& g2; // ions

   group_t& product1;
   group_t& product2;
   CollisionSpec params;

   particle_vec g1_products; // electrons
   particle_vec g2_products; // ions

   std::array<std::mt19937_64, nThreads> generator;

   bool has_coulomb;
   bool has_ionization;

   Collisions(group_t& g1_, group_t& g2_, group_t& prod1_, group_t& prod2_, const auto& params_)
   : g1(g1_),
     g2(g2_),
     product1(prod1_),
     product2(prod2_),
     params(params_),
     has_coulomb(std::ranges::contains(params_.channels, "coulomb")),
     has_ionization(std::ranges::contains(params_.channels, "ionization"))
   {
      for (std::size_t i = 0; i < nThreads; i++) {
         generator[i] = init_mt_64();
      }
   }

   static std::mt19937_64 init_mt_64() {
      std::array<int, 624> seed_data{};
      std::random_device r;
      std::generate_n(seed_data.data(), seed_data.size(), std::ref(r));
      std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
      return std::mt19937_64(seq);
   }



   void advance(const auto step) requires (coll_enabled) {
      static constexpr auto cell_vol_inv = 1.0 / (dx * dy * dz);
      static constexpr auto dt_vol = dt * cell_vol_inv;
      static constexpr auto q_eps = constants::q_e * constants::eps0;
      static constexpr auto pi43 = 4.0 * constants::pi / 3.0;
      static constexpr auto twothirds = 2.0 / 3.0;

      const auto pi34_cuberoot = std::pow(pi43, 1.0 / 3.0);
      constexpr auto four_pi_eps_c2_inv = 1.0 / (4.0 * constants::pi * constants::eps0 * constants::c_sqr);

      // todo: these can probably be moved out of this function
      const auto q1q2 = g1.charge * g2.charge;
      const auto m1m2 = g1.mass * g2.mass;
      const auto m1_over_m2 = g1.mass / g2.mass;
      const auto m1c2 = g1.charge * constants::c_sqr;
      const auto m2c2 = g2.charge * constants::c_sqr;

      const auto coef1 = four_pi_eps_c2_inv * math::SQR(q1q2) / (m1m2 * constants::eps0);
      const auto coef2 = four_pi_eps_c2_inv * std::abs(q1q2);

      if (step % params.step_interval != 0) { return; }

      g1.update_cell_map();
      g2.update_cell_map();

      g1_products.clear();
      g2_products.clear();

      for (auto& [z_code, cell1] : g1.cell_map) {
         if (!g2.cell_map.contains(z_code)) {
            continue;
         }

         const auto& cell2 = g2.cell_map.at(z_code);
         const auto np1 = cell1.size();
         const auto np2 = cell2.size();
         const auto np1_lt_np2 = np1 < np2;

         const auto n_partners = params.self_scatter ? np1 - 1 + np1 % 2 : std::max(np1, np2);
         const auto scatter_coef = dt_vol * static_cast<double>(n_partners);

         // todo: move vector init outside of loop
         std::uniform_real_distribution rng(0.0, 1.0);
         std::vector<std::size_t> pids1(n_partners);
         std::vector<std::size_t> pids2(n_partners);
         std::vector<double> nDups(n_partners);

         CoulombData cell_data{};
         if (has_coulomb) {
            // todo: Cell data here is really only used for Coulomb
            auto calcDensityTemp = [](const auto& c, const auto& mc2) -> std::tuple<double, double> {
               // Calculates cell density and temperature in eV
               auto ttl_weight = 0.0;
               auto KE = 0.0;
               for (auto& p : c) {
                  ttl_weight += p.weight;
                  KE += (p.gamma - 1.0) * p.weight * mc2;
               }
               return {ttl_weight * cell_vol_inv, twothirds * KE / (ttl_weight * constants::q_e)};
            };

            const auto& rhoT1 = calcDensityTemp(cell1, m1c2);
            const auto& [density1, temp1] = rhoT1;
            const auto& [density2, temp2] = params.self_scatter ? rhoT1 : calcDensityTemp(cell2, m2c2);

            const auto rmin2 = std::pow(pi43 * std::max(density1, density2), -twothirds);
            const auto lD2_1 = q_eps * temp1 / (density1 * math::SQR(g1.charge));
            const auto lD2_2 = params.self_scatter ? 0.0 : q_eps * temp2 / (density2 * math::SQR(g2.charge));
            const auto bmax2 = rmin2 + lD2_1 + lD2_2;

            const auto d1 = m1_over_m2 * std::pow(density1, twothirds);
            const auto d2 = std::pow(density2, twothirds);
            const auto scatter_lowT = scatter_coef * pi34_cuberoot * (m1_over_m2 + 1.0) / std::max(d1, d2);

            cell_data = {coef1, coef2, bmax2, scatter_lowT};
         }

         // Create vector of random integers in range [0, Np1) and [0, Np2)
         for (std::size_t i = 0; i < n_partners; ++i) {
            pids1[i] = i % np1;
            pids2[i] = i % np2;
         }

         const auto grp_dups = static_cast<double>(np1_lt_np2 ? np2 / np1 : np1 / np2);
         const auto      rem =    static_cast<int>(np1_lt_np2 ? np2 % np1 : np1 % np2);
         std::ranges::fill(nDups, grp_dups);
         if (rem != 0) {
            std::ranges::fill(nDups.end() - rem, nDups.end(), grp_dups);
         }

         std::ranges::shuffle(pids1, generator[0]);
         std::ranges::shuffle(pids2, generator[0]);

         #pragma omp parallel for num_threads(nThreads)
         for (std::size_t i = 0; i < n_partners; ++i) {
            const auto tid = static_cast<std::size_t>(omp_get_thread_num());
            const auto pid1 = pids1[i];
            const auto pid2 = pids2[i];

            const auto dup1 =  np1_lt_np2 ? nDups[pid1] : 1.0;
            const auto dup2 = !np1_lt_np2 ? nDups[pid2] : 1.0;
            const auto weight1 = cell1[pid1].weight / dup1;
            const auto weight2 = cell1[pid2].weight / dup2;
            const auto max_weight = std::max(weight1, weight2);

            ParticlePairData pair_data {
               cell1[pid1],
               cell2[pid2],
               g1.mass,
               g2.mass,
               weight1,
               weight2,
               max_weight,
               scatter_coef,
               {rng(generator[tid]), rng(generator[tid]), rng(generator[tid]), rng(generator[tid])}
            };

            coulombCollision(has_coulomb, pair_data, params, cell_data);

            // todo: add randomly selecting channel when more channels are added.
            ionizationCollision(
               has_ionization,
               pair_data,
               params,
               product1,
               product2,
               nDups[pids1[i]] // todo: should this be pid1[i] or pid2[i] or min(pid1, pid2)?
            );
         } // end for(npairs)
      } // end for(z_code, cell1)

      if (!g1_products.empty()) {
         // G1 should always be electrons
         product1.particles.insert(product1.particles.end(), g1_products.begin(), g1_products.end());
         g1.cell_map_updated = false;
         g1.is_sorted = false;
         g1.sort_particles();
      }

      if (!g2_products.empty()) {
         product2.particles.insert(product2.particles.end(), g2_products.begin(), g2_products.end());
         g2.cell_map_updated = false;
         g2.is_sorted = false;
         g2.sort_particles();
      }
   } // end update()

   static void advance(const auto) requires (!coll_enabled) {}

}; // end struct Collisions

} // end namespace tf::collisions

#endif //COLLISIONS_HPP
