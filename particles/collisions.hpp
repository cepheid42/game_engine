#ifndef COLLISIONS_HPP
#define COLLISIONS_HPP

#include "binary_channels.hpp"
#include "constants.hpp"
#include "interpolation.hpp"
#include "math_utils.hpp"
#include "particles.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <omp.h>
#include <random>
#include <ranges>
#include <tuple>
#include <vector>

namespace tf::collisions
{

struct Collisions {
   using group_t = particles::ParticleGroup;
   using particle_vec = std::vector<particles::Particle>;

   struct Products {
      using group_t = particles::ParticleGroup;
      group_t* product1{nullptr};
      group_t* product2{nullptr};
   };

   struct Buffers {
      particle_vec g1_products{};
      particle_vec g2_products{};
   };

   group_t& g1;
   group_t& g2;
   std::map<std::string, Products> products{};
   std::map<std::string, Buffers> buffers{};

   CollisionSpec specs;

   static constexpr auto NRNG = 2000;
   std::mt19937_64 generator;
   std::uniform_real_distribution<> rng;
   std::vector<double> rngs;

   bool has_coulomb;
   enum struct ChannelType { None, Ionization, Fusion, Bremmstrahlung };
   std::vector<ChannelType> channels{};

   interp::Table ionization_cs{};
   interp::Table fusion_cs{};
   interp::MultiTable brem_cs{};

   Collisions(const auto& params_, auto& group_map)
   : g1(group_map.at(std::string{params_.group1})),
     g2(group_map.at(std::string{params_.group2})),
     specs(params_),
     generator(init_mt_64()),
     rng(0.0, 1.0),
     rngs(NRNG),
     has_coulomb(std::ranges::contains(params_.channels, "coulomb"))
   {
      if (std::ranges::contains(params_.channels, "ionization")) {
         channels.push_back(ChannelType::Ionization);
         products.emplace(
               "ionization",
               Products{&group_map.at(std::string{params_.ionization.product1}),
                        &group_map.at(std::string{params_.ionization.product2})}
         );
         buffers.emplace("ionization", Buffers{});

         if (not specs.ionization.cross_section_file.empty() and specs.ionization.constant_cross_section == 0.0) {
            ionization_cs = interp::Table(std::string{specs.ionization.cross_section_file}, 1.0, 1.0); // eV and m^2
         }
      }

      if (std::ranges::contains(params_.channels, "fusion")) {
         channels.push_back(ChannelType::Fusion);
         products.emplace(
            "fusion",
            Products{&group_map.at(std::string(params_.fusion.product1)),
                     &group_map.at(std::string(params_.fusion.product2))}
         );
         buffers.emplace("fusion", Buffers{});

         if (not specs.fusion.cross_section_file.empty() and specs.fusion.constant_cross_section == 0.0) {
            fusion_cs = interp::Table(std::string{specs.fusion.cross_section_file}, 1.0, 1.0);
         }
      }

      if (std::ranges::contains(params_.channels, "radiation")) {
         channels.push_back(ChannelType::Bremmstrahlung);
         products.emplace(
            "radiation",
            Products{&group_map.at(std::string(params_.radiation.product1))}
         );
         buffers.emplace("radiation", Buffers{});

         if (not specs.radiation.cross_section_file.empty()) {
            brem_cs = interp::MultiTable(std::string{specs.radiation.cross_section_file});
         }
      }

      if (channels.empty()) { channels.push_back(ChannelType::None); }
   }

   static std::mt19937_64 init_mt_64() {
      std::array<int, 624> seed_data{};
      std::random_device r;
      std::generate_n(seed_data.data(), seed_data.size(), std::ref(r));
      std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
      return std::mt19937_64(seq);
   }

   void applyCollision(const ChannelType type, const ParticlePairData& pair_data) {
      switch (type)
      {
         case ChannelType::Ionization:
            {
               ionizationCollision(
                  pair_data,
                  specs,
                  buffers["ionization"],
                  ionization_cs
               );
               break;
            }
         case ChannelType::Fusion:
            {
               fusionCollision(
                  pair_data,
                  specs,
                  buffers["fusion"],
                  fusion_cs,
                  products["fusion"].product1->mass,
                  products["fusion"].product2->mass
               );
               break;
            }
         case ChannelType::Bremmstrahlung:
            {
               bremsstrahlungCollision(
                  pair_data,
                  specs,
                  buffers["radiation"],
                  brem_cs
               );
            }
         default:
            break;
      }
   }

   CoulombData getCoulombData(const auto& cell1, const auto& cell2, const auto scatter_coef) {
      static constexpr auto cell_vol_inv = 1.0 / (dx * dy * dz);
      static constexpr auto q_eps = constants::q_e * constants::eps0;
      static constexpr auto four_pi_eps_c2_inv = 1.0 / (4.0 * constants::pi * constants::eps0 * constants::c_sqr);
      static const auto pi34_cuberoot = std::pow((4.0 / 3.0) * constants::pi, 1.0 / 3.0);

      if (not has_coulomb) { return {}; } // No need to calculate all this if no coulomb collisions

      const auto q1q2 = g1.charge * g2.charge;
      const auto m1_over_m2 = g1.mass / g2.mass;
      const auto coef1 = four_pi_eps_c2_inv * math::SQR(q1q2) / (g1.mass * g2.mass * constants::eps0);
      const auto coef2 = four_pi_eps_c2_inv * std::abs(q1q2);

      auto calcDensityTemp = [](const auto& c, const auto& mc2) -> std::tuple<double, double> {
         // Calculates cell density and temperature in eV
         auto ttl_weight = 0.0;
         auto KE = 0.0;
         for (auto& p : c) {
            ttl_weight += p.weight;
            KE += (p.gamma - 1.0) * p.weight * mc2;
         }
         return {ttl_weight * cell_vol_inv, (2.0 / 3.0) * KE / (ttl_weight * constants::q_e)};
      };

      const auto& rhoT1 = calcDensityTemp(cell1, g1.mass * constants::c_sqr);
      const auto& [density1, temp1] = rhoT1;
      const auto& [density2, temp2] = specs.self_scatter ? rhoT1 : calcDensityTemp(cell2, g2.mass * constants::c_sqr);

      const auto rmin2 = std::pow((4.0 / 3.0) * constants::pi * std::max(density1, density2), -(2.0 / 3.0));
      const auto lD2_1 = q_eps * temp1 / (density1 * math::SQR(g1.charge));
      const auto lD2_2 = specs.self_scatter ? 0.0 : q_eps * temp2 / (density2 * math::SQR(g2.charge));
      const auto bmax2 = rmin2 + lD2_1 + lD2_2;

      const auto d1 = m1_over_m2 * std::pow(density1, (2.0 / 3.0));
      const auto d2 = std::pow(density2, (2.0 / 3.0));
      const auto scatter_lowT = scatter_coef * pi34_cuberoot * (m1_over_m2 + 1.0) / std::max(d1, d2);

      return {coef1, coef2, bmax2, scatter_lowT};
   }

   void advance(const auto step) requires (coll_enabled) {
      if (step % specs.step_interval != 0) { return; }

      const auto reduced_mass = (g1.mass * g2.mass) / (g1.mass + g2.mass);

      g1.cell_map_updated = false;
      g1.is_sorted = false;
      g1.update_cell_map();

      g2.cell_map_updated = false;
      g2.is_sorted = false;
      g2.update_cell_map();

      for (auto& [g1_products, g2_products] : buffers | std::views::values) {
         g1_products.clear();
         g2_products.clear();
      }

      std::ranges::generate(rngs, [&]{ return rng(generator); });
      
      std::vector<std::size_t> cell_ids(g1.cell_map.size());
      std::transform(g1.cell_map.begin(), g1.cell_map.end(), cell_ids.begin(), [](auto& kv) { return kv.first; });

      #pragma omp parallel for num_threads(nThreads)
      for (auto j = 0zu; j < cell_ids.size(); j++) {
         const auto tid = omp_get_thread_num();
         const auto z_code = cell_ids[j];

         if (not (g1.cell_map.contains(z_code) and g2.cell_map.contains(z_code))) {
            continue;
         }
      
         const auto& cell1 = g1.cell_map.at(z_code);
         const auto& cell2 = g2.cell_map.at(z_code);
      
         const auto np1 = cell1.size();
         const auto np2 = cell2.size();
         assert(np1 != 0 and np2 != 0);
      
         const auto np1_lt_np2 = np1 < np2;
      
         const auto n_partners = specs.self_scatter ? np1 - 1 + np1 % 2 : std::max(np1, np2);
         const auto n_pairs = specs.self_scatter ? n_partners / 2 : n_partners;
         const auto scatter_coef = static_cast<double>(n_partners) * dt / (dx * dy * dz);
      
         std::vector<std::size_t> pids1(n_partners);
         std::vector<std::size_t> pids2(n_partners);
         std::vector<double> nDups(std::min(np1, np2));

         const auto cell_data = getCoulombData(cell1, cell2, scatter_coef);

         // Create vector of random integers in range [0, Np1) and [0, Np2)
         for (auto i = 0zu; i < n_partners; ++i) {
            pids1[i] = i % np1;
            pids2[i] = i % np2;
         }
      
         const auto grp_dups = static_cast<double>(np1_lt_np2 ? np2 / np1 : np1 / np2);
         const auto      rem =    static_cast<int>(np1_lt_np2 ? np2 % np1 : np1 % np2);
         std::ranges::fill(nDups, grp_dups);
         if (rem != 0) {
            std::ranges::fill(nDups.begin(), nDups.begin() + rem, grp_dups + 1);
         }
      
         std::ranges::shuffle(pids1, generator);
         // std::ranges::shuffle(pids2, generator);

         for (auto i = 0zu; i < n_pairs; ++i) {
            const auto pid1 = pids1[i];
            const auto pid2 = pids2[i];

            auto& p1 = cell1[pid1];
            auto& p2 = cell2[pid2];
            if (p1.is_disabled() or p2.is_disabled() or (specs.self_scatter and pid1 == pid2)) { continue; }
      
            const auto dups = std::max(
                np1_lt_np2 ? nDups[pid1] : 1.0,
               !np1_lt_np2 ? nDups[pid2] : 1.0
            );
            const auto weight1 = p1.weight / dups;
            const auto weight2 = p2.weight / dups;

            const auto idx = (static_cast<int>(i) + tid) % (NRNG - 6);

            ParticlePairData pair_data {
               p1, p2, reduced_mass, g1.mass, g2.mass, dups, scatter_coef,
               std::max(weight1, weight2),
               std::span{rngs.begin() + idx, rngs.begin() + idx + 5}
            };

            if (has_coulomb) { coulombCollision(pair_data, specs, cell_data); }

            const auto channel_idx = static_cast<std::size_t>(rngs[0] * static_cast<double>(channels.size()));
            assert(channel_idx <= channels.size() - 1);
            applyCollision(channels[channel_idx], pair_data);
         } // end for(npairs)
      } // end for(z_code, cell1)

      for (auto& [k, v] : buffers) {
         // todo: check this logic to make sure its valid for all collision types (not just ionization)
         if (not v.g1_products.empty()) {
            auto& product1 = *(products[k].product1);
            product1.particles.insert(product1.particles.end(), v.g1_products.begin(), v.g1_products.end());
            product1.cell_map_updated = false;
            product1.is_sorted = false;
         }

         if (not v.g2_products.empty()) {
            auto& product2 = *(products[k].product2);
            product2.particles.insert(product2.particles.end(), v.g2_products.begin(), v.g2_products.end());
            product2.cell_map_updated = false;
            product2.is_sorted = false;
         }
      } // end for(buffers)
   } // end update()

   static void advance(const auto) requires (not coll_enabled) {}

}; // end struct Collisions

} // end namespace tf::collisions

#endif //COLLISIONS_HPP
