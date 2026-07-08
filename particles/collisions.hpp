#ifndef COLLISIONS_HPP
#define COLLISIONS_HPP

#include "binary_channels.hpp"
#include "constants.hpp"
#include "interpolation.hpp"
#include "math_utils.hpp"
#include "particles.hpp"
#include "timers.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <map>
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

   enum struct ChannelType { None, Ionization, FusionA, FusionB, Bremsstrahlung, BremsstrahlungTFD };

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
   std::map<ChannelType, Products> products{};
   std::map<ChannelType, Buffers> buffers{};
   std::map<ChannelType, interp::Table> tables{};

   CollisionSpec specs;

   static constexpr auto NRNG = 2000;
   std::mt19937_64 generator;
   std::uniform_real_distribution<> rng;
   std::vector<double> rngs;

   bool has_coulomb;
   bool has_tfd_radiation;
   std::vector<ChannelType> channels{};

   interp::BremTFDTable bremtfd_cs{};
   interp::BremTable brem_cs{};

   utilities::Timer step_timer{};

   Collisions(const auto& params_, auto& group_map)
   : g1(group_map.at(std::string{params_.group1})),
     g2(group_map.at(std::string{params_.group2})),
     specs(params_),
     generator(init_mt_64()),
     rng(0.0, 1.0),
     rngs(NRNG),
     has_coulomb(std::ranges::contains(params_.channels, "coulomb")),
     has_tfd_radiation(std::ranges::contains(params_.channels, "radiation") and params_.radiation.use_TFD)
   {
      if (std::ranges::contains(params_.channels, "ionization")) {
         channels.push_back(ChannelType::Ionization);
         products.emplace(
               ChannelType::Ionization,
               Products{&group_map.at(std::string{params_.ionization.product1}),
                        &group_map.at(std::string{params_.ionization.product2})}
         );
         buffers.emplace(ChannelType::Ionization, Buffers{});

         if (not specs.ionization.cross_section_file.empty() and specs.ionization.constant_cross_section == 0.0) {
            tables.emplace(ChannelType::Ionization, interp::Table(std::string{specs.ionization.cross_section_file}, 1.0, 1.0)); // eV and m^2
         }
      }

      if (std::ranges::contains(params_.channels, "fusion")) {
         auto channel = 0;
         for (auto& spec : params_.fusion) {
            if (spec.product1.empty() and spec.product2.empty()) {
               continue;
            }

            const auto type = (channel == 0) ? ChannelType::FusionA : ChannelType::FusionB;
            channels.push_back(type);
            products.emplace(
               type,
               Products{&group_map.at(std::string(spec.product1)),
                        &group_map.at(std::string(spec.product2))}
            );

            buffers.emplace(type, Buffers{});

            if (not spec.cross_section_file.empty() and spec.constant_cross_section == 0.0) {
               tables.emplace(type, interp::Table(std::string{spec.cross_section_file}, 1.0, 1.0));
            }
            channel++;
         }
      }

      if (std::ranges::contains(params_.channels, "radiation")) {
         const auto type = params_.radiation.use_TFD ? ChannelType::BremsstrahlungTFD : ChannelType::Bremsstrahlung;
         channels.push_back(type);
         products.emplace(
            type,
            Products{&group_map.at(std::string(params_.radiation.product1))}
         );
         buffers.emplace(type, Buffers{});

         if (params_.radiation.use_TFD) {
            bremtfd_cs = interp::BremTFDTable(
               g2.atomic_number,
               g2.charge,
               specs.radiation.min_energy,
               specs.radiation.max_energy
            );
         }


         if (not params_.radiation.use_TFD and not specs.radiation.cross_section_file.empty()) {
            brem_cs = interp::BremTable(std::string{specs.radiation.cross_section_file});
         }
      }

      if (channels.empty()) { channels.push_back(ChannelType::None); }
   }

   ~Collisions() {
      std::println("{}-{} Collisions: {}", g1.name, g2.name, std::chrono::hh_mm_ss(step_timer.elapsed));
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
                  specs.ionization,
                  specs.probability_search_area,
                  buffers[ChannelType::Ionization],
                  tables[ChannelType::Ionization]
               );
               break;
            }
         case ChannelType::FusionA:
            {
               fusionCollision(
                  pair_data,
                  specs.fusion[0],
                  specs.probability_search_area,
                  buffers[ChannelType::FusionA],
                  tables[ChannelType::FusionA],
                  products[ChannelType::FusionA].product1->mass,
                  products[ChannelType::FusionA].product2->mass
               );
               break;
            }
         case ChannelType::FusionB:
            {
               fusionCollision(
                  pair_data,
                  specs.fusion[1],
                  specs.probability_search_area,
                  buffers[ChannelType::FusionB],
                  tables[ChannelType::FusionB],
                  products[ChannelType::FusionB].product1->mass,
                  products[ChannelType::FusionB].product2->mass
               );
               break;
            }
         case ChannelType::Bremsstrahlung:
            {
               bremsstrahlungCollision(
                  pair_data,
                  specs.radiation,
                  buffers[ChannelType::Bremsstrahlung],
                  brem_cs
               );
               break;
            }
         case ChannelType::BremsstrahlungTFD:
            {
               bremsstrahlungCollision(
                  pair_data,
                  specs.radiation,
                  buffers[ChannelType::BremsstrahlungTFD],
                  bremtfd_cs
               );
               break;
            }
         default:
            break;
      }
   }

   CellData getCellData(const auto& cell1, const auto& cell2, const auto cid, const auto scatter_coef) {
      static constexpr auto cell_vol_inv = 1.0 / (dx * dy * dz);
      static constexpr auto q_eps = constants::q_e * constants::eps0;
      static constexpr auto four_pi_eps_c2_inv = 1.0 / (4.0 * constants::pi * constants::eps0 * constants::c_sqr);
      static const auto pi34_cuberoot = std::pow((4.0 / 3.0) * constants::pi, 1.0 / 3.0);

      if (not (has_coulomb or has_tfd_radiation)) { return {}; } // No need to calculate all this if no coulomb collisions

      const auto q1q2 = g1.charge * g2.charge;
      const auto m1_over_m2 = g1.mass / g2.mass;
      const auto coef1 = four_pi_eps_c2_inv * math::SQR(q1q2) / (g1.mass * g2.mass * constants::eps0);
      const auto coef2 = four_pi_eps_c2_inv * std::abs(q1q2);

      auto calcDensityTemp = [](const auto& c, const auto& mc2) -> std::tuple<double, double> {
         // Calculates cell density and temperature in eV.
         // The temperature moment is also used to build a cell-averaged Coulomb log.
         auto ttl_weight = 0.0;
         auto KE = 0.0;
         for (auto& p : c) {
            ttl_weight += p.weight;
            KE += (p.gamma - 1.0) * p.weight * mc2;
         }
         if (ttl_weight <= 0.0) { return {0.0, 0.0}; }
         return {ttl_weight * cell_vol_inv, (2.0 / 3.0) * KE / (ttl_weight * constants::q_e)};
      };

      const auto& rhoT1 = calcDensityTemp(cell1, g1.mass * constants::c_sqr);
      const auto& [density1, temp1] = rhoT1;
      const auto& [density2, temp2] = specs.self_scatter ? rhoT1 : calcDensityTemp(cell2, g2.mass * constants::c_sqr);

      if (density1 <= 0.0 or density2 <= 0.0) { return {}; }

      const auto rmin2 = std::pow((4.0 / 3.0) * constants::pi * std::max(density1, density2), -(2.0 / 3.0));
      const auto lD2_1 = q_eps * std::max(temp1, 1.0e-30) / (density1 * math::SQR(g1.charge));
      const auto lD2_2 = specs.self_scatter ? 0.0 : q_eps * std::max(temp2, 1.0e-30) / (density2 * math::SQR(g2.charge));
      const auto bmax2 = rmin2 + lD2_1 + lD2_2;

      const auto d1 = m1_over_m2 * std::pow(density1, (2.0 / 3.0));
      const auto d2 = std::pow(density2, (2.0 / 3.0));
      const auto scatter_lowT = scatter_coef * pi34_cuberoot * (m1_over_m2 + 1.0) / std::max(d1, d2);

      auto coulomb_log_cell = specs.coulomb.coulomb_log;
      if (has_coulomb and coulomb_log_cell <= 0.0) {
         // Cell-averaged Coulomb log: this avoids recomputing a logarithm and de Broglie
         // length for every Monte Carlo pair.  The reference relative velocity is thermal
         // and is intentionally floored so cold cells remain finite.
         const auto reduced_mass = (g1.mass * g2.mass) / (g1.mass + g2.mass);
         const auto vrel2_ref = std::max(
            1.0,
            3.0 * std::max(temp1, 1.0e-30) * constants::q_e / g1.mass +
            (specs.self_scatter ? 3.0 * std::max(temp1, 1.0e-30) * constants::q_e / g1.mass
                                : 3.0 * std::max(temp2, 1.0e-30) * constants::q_e / g2.mass)
         );
         const auto pstar_ref = std::max(1.0e-300, reduced_mass * std::sqrt(vrel2_ref));
         const auto l_deBroglie = 0.5 * constants::h / pstar_ref;
         const auto gamma_coef1_nr = 1.0 / (g1.mass + g2.mass);
         const auto gamma_coef2_nr = 1.0 + constants::c_sqr * g1.mass * g2.mass / math::SQR(pstar_ref);
         const auto b0 = gamma_coef1_nr * gamma_coef2_nr * coef2;
         const auto bmin2 = std::max(math::SQR(b0), math::SQR(l_deBroglie));
         coulomb_log_cell = std::max(2.0, 0.5 * std::log(1.0 + bmax2 / bmin2));
      }

      const auto l_debye = std::sqrt(lD2_1 + lD2_2);
      if (has_tfd_radiation and l_debye != 0.0) {
         bremtfd_cs.updateCDF(cid, std::max(l_debye, std::sqrt(rmin2)));
      }

      return {coef1, coef2, bmax2, scatter_lowT, coulomb_log_cell};
   }

   void advance(const auto step) requires (coll_enabled) {
      step_timer.start_timer();
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

      std::vector<std::mt19937_64> thread_generators(nThreads);
      for (auto& gen : thread_generators) { gen.seed(generator()); }

      std::vector<std::size_t> cell_ids(g1.cell_map.size());
      std::transform(g1.cell_map.begin(), g1.cell_map.end(), cell_ids.begin(), [](auto& kv) { return kv.first; });

      const auto has_only_coulomb = has_coulomb and channels.size() == 1 and channels[0] == ChannelType::None;

      #pragma omp parallel num_threads(nThreads)
      {
         const auto tid = omp_get_thread_num();
         auto& local_generator = thread_generators[tid];
         std::uniform_real_distribution<> local_rng(0.0, 1.0);

         std::vector<double> local_rngs(NRNG);
         std::vector<std::size_t> pids1(100);
         std::vector<std::size_t> pids2(100);
         std::vector<double> nDups(25);

         #pragma omp for
         for (auto j = 0zu; j < cell_ids.size(); j++) {
            const auto z_code = cell_ids[j];

            if (not (g1.cell_map.contains(z_code) and g2.cell_map.contains(z_code))) {
               continue;
            }

            const auto& cell1 = g1.cell_map.at(z_code);
            const auto& cell2 = g2.cell_map.at(z_code);

            const auto np1 = cell1.size();
            const auto np2 = cell2.size();
            if (np1 == 0 or np2 == 0) { continue; }
            assert(np1 != 0 and np2 != 0);

            const auto np1_lt_np2 = np1 < np2;

            const auto n_partners = specs.self_scatter ? np1 - 1 + np1 % 2 : std::max(np1, np2);
            const auto n_pairs = specs.self_scatter ? n_partners / 2 : n_partners;
            const auto scatter_coef = static_cast<double>(n_partners) * dt / (dx * dy * dz);

            pids1.resize(n_partners);
            pids2.resize(n_partners);
            nDups.resize(std::min(np1, np2));

            const auto cell_data = getCellData(cell1, cell2, z_code, scatter_coef);

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

            // Use a thread-local generator here.  The old shared generator inside the
            // OpenMP loop was both slower and not thread safe.
            std::ranges::shuffle(pids1, local_generator);
            // std::ranges::shuffle(pids2, local_generator);

            std::ranges::generate(local_rngs, [&]{ return local_rng(local_generator); });

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

               const auto idx = static_cast<int>(i) % (NRNG - 6);
               auto pair_rng = std::span{local_rngs.begin() + idx, local_rngs.begin() + idx + 5};

               ParticlePairData pair_data {
                  p1, p2, reduced_mass, g1.mass, g2.mass, dups, scatter_coef,
                  std::max(weight1, weight2), z_code,
                  pair_rng
               };

               if (has_coulomb) { coulombCollision(pair_data, specs, cell_data); }

               // In the common Coulomb-only case, skip the channel dispatch and map lookups.
               if (has_only_coulomb) { continue; }

               const auto channel_idx = static_cast<std::size_t>(pair_rng[3] * static_cast<double>(channels.size()));
               assert(channel_idx <= channels.size() - 1);
               applyCollision(channels[channel_idx], pair_data);
            } // end for(npairs)
         } // end for(z_code, cell1)
      } // end omp parallel

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
      step_timer.stop_timer();
   } // end update()

   static void advance(const auto) requires (not coll_enabled) {}

}; // end struct Collisions
} // end namespace tf::collisions

#endif //COLLISIONS_HPP
