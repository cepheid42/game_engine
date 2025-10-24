#ifndef COLLISIONS_HPP
#define COLLISIONS_HPP

#include "particles.hpp"
#include "constants.hpp"
#include "math_utils.hpp"
#include "rng_utils.h"

// #include <omp.h>

#include <algorithm>
#include <map>
#include <numeric>
#include <span>
#include <print>

namespace tf::collisions
{
// enum class EmissionType { CoulombSmallAngle, Isotropic, FusionPolyFit };
enum class CollisionType {
   Coulomb,
   // Fusion,
   // Ionization
   // Neutral,
   // Excitation,
   // Compton,
   // PhotonAbsorption,
   // Brehmsstrahlung (with or without TFD)
};

inline constexpr std::array Park2022_ScatterCoefs = {
   std::array{0.283,  0.667,  0.0307, 16.971,  6.59, 29.02, 0.0258, 0.295, 0.00328, 1.794}, // He
   std::array{0.0502, 0.9190, 0.0799,  0.9521, 3.29, 10.91, 0.008,  0.936, 0.0177,  1.251}, // H, D, T
   std::array{0.409,  0.5737, 0.2665,  3.912,  2.81,  4.33, 0.006,  0.0,   0.0,     0.729}  // H2, D2, T2
};

// template<>
struct Collisions {
   using group_t = particles::ParticleGroup;

   group_t& g1;
   group_t& g2;

   // std::vector<std::size_t> pids1{};
   // std::vector<std::size_t> pids2{};
   // std::vector<double> nDups{};

   double coulombLog;
   double rate_mult;
   int step_interval;
   bool self_scatter;

   std::array<std::mt19937_64, nThreads> generator;

   Collisions(group_t& g1_, group_t& g2_, const auto clog_, const auto mult_, const auto step_, const auto self_scatter_)
   : g1(g1_), g2(g2_), coulombLog(clog_), rate_mult(mult_), step_interval(step_), self_scatter(self_scatter_)
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

   auto calculate_P1final(const auto& p1_star, const auto& vcm, const auto vcm_len2, const auto gamma_cm, const auto rand0, const auto getScatteringAngle) {
      const auto p1_star_len = p1_star.length();

      const auto [cos_theta, sin_theta] = getScatteringAngle();
      const auto phi = 2.0 * constants::pi<double> * rand0;
      const auto dv0 = sin_theta * std::cos(phi);
      const auto dv1 = sin_theta * std::sin(phi);

      const auto p1_star_perp = std::sqrt(math::SQR(p1_star[0]) + math::SQR(p1_star[1]));
      const vec3 p1f_star {
         dv0 * (p1_star[0] * p1_star[2] / p1_star_perp) - dv1 * (p1_star[1] * p1_star_len / p1_star_perp) + cos_theta * p1_star[0],
         dv0 * (p1_star[1] * p1_star[2] / p1_star_perp) + dv1 * (p1_star[0] * p1_star_len / p1_star_perp) + cos_theta * p1_star[1],
         -dv0 * p1_star_perp + cos_theta * p1_star[2]
      };

      const auto pcoef_lab = (gamma_cm - 1.0) * dot(vcm, p1f_star) / vcm_len2;

      return {p1f_star, pcoef_lab};
   }

   auto applyScatter(auto& particle1, auto& particle2, const auto& vcm, const auto& p1f_star, const auto pcoef_lab, const auto gamma1_star, const auto gamma2_star, const auto gamma_cm, const auto condition1, const auto condition2) {
      if (condition1) {
         const auto p1f = p1f_star + (pcoef_lab + g1.mass * gamma1_star * gamma_cm) * vcm;
         const auto gamma1f = std::sqrt(1.0 + p1f.length_squared() * constants::over_c_sqr<double> / math::SQR(g1.mass));
         particle1.velocity = p1f / (g1.mass * gamma1f);
         particle1.gamma = gamma1f;
      }

      if (condition2) {
         const auto p2f = -p1f_star + (-pcoef_lab + g2.mass * gamma2_star * gamma_cm) * vcm;
         const auto gamma2f = std::sqrt(1.0 + p2f.length_squared() * constants::over_c_sqr<double> / math::SQR(g2.mass));
         particle2.velocity = p2f / (g2.mass * gamma2f);
         particle2.gamma = gamma2f;
      }
   }

   struct COMData {
      double mg1;
      double mg2;
      double gamma_cm;
      double gamma1_cm;
      double gamma2_cm;
      double vcm_dot_v1;
      double vcm_dot_v2;
      vec3<double> vcm;
      vec3<double> p1;
      vec3<double> p2;
      vec3<double> p1_star;
   };

   COMData getCOMData(const auto& particle1, const auto& particle2, const auto m1, const auto m2) {
      const auto mg1 = m1 * particle1.gamma;
      const auto mg2 = m2 * particle2.gamma;

      const auto p1 = mg1 * particle1.velocity;
      const auto p2 = mg2 * particle2.velocity;

      const auto vcm = (p1 + p2) / (mg1 + mg2);
      const auto vcm_len2 = vcm.length_squared();
      const auto gamma_cm = particles::calculateGammaV(vcm);

      const auto vcm_dot_v1 = dot(vcm, particle1.velocity);
      const auto vcm_dot_v2 = dot(vcm, particle2.velocity);

      const auto p1_star = p1 + vcm * mg1 * (vcm_dot_v1 / vcm_len2 * (gamma_cm - 1.0) - gamma_cm);

      const auto gamma1_cm = gamma_cm * (1.0 - vcm_dot_v1 * constants::over_c_sqr<double>);
      const auto gamma2_cm = gamma_cm * (1.0 - vcm_dot_v2 * constants::over_c_sqr<double>);

      return {
         .mg1 = mg1,
         .mg2 = mg2,
         .gamma_cm = gamma_cm,
         .gamma1_cm = gamma1_cm,
         .gamma2_cm = gamma2_cm,
         .vcm_dot_v1 = vcm_dot_v1,
         .vcm_dot_v2 = vcm_dot_v2,
         .vcm = vcm,
         .p1 = p1,
         .p2 = p2,
         .p1_star = p1_star,
      };
   }

   void coulombCollision(auto& particle1, auto& particle2,
      const auto weight1, const auto weight2, const auto max_weight,
      const auto rand0, const auto rand1, const auto rand2,
      const auto coef1, const auto coef2, const auto bmax2,
      const auto scatter_coef, const auto scatter_lowT)
   {
      static constexpr auto h_half = 0.5 * constants::h<double>;

      const auto scatter_p1 = rand0 * weight1 < weight2;
      const auto scatter_p2 = rand0 * weight2 < weight1;
      const auto dv_length = (particle1.velocity - particle2.velocity).length();

      if (dv_length == 0.0 or !(scatter_p1 or scatter_p2)) {
         return;
      } // parallel particles or self particles can cause divide by zero errors

      // const auto mg1 = g1.mass * particle1.gamma;
      // const auto mg2 = g2.mass * particle2.gamma;
      //
      // const auto p1 = mg1 * particle1.velocity;
      // const auto p2 = mg2 * particle2.velocity;
      //
      // const auto vcm = (p1 + p2) / (mg1 + mg2);
      // const auto vcm_len2 = vcm.length_squared();
      // const auto gamma_cm = particles::calculateGammaV(vcm);
      //
      // const auto vcm_dot_v1 = dot(vcm, particle1.velocity);
      // const auto vcm_dot_v2 = dot(vcm, particle2.velocity);
      //
      // const auto p1_star = p1 + vcm * mg1 * (vcm_dot_v1 / vcm_len2 * (gamma_cm - 1.0) - gamma_cm);
      //
      // const auto gamma1_cm = gamma_cm * (1.0 - vcm_dot_v1 * constants::over_c_sqr<double>);
      // const auto gamma2_cm = gamma_cm * (1.0 - vcm_dot_v2 * constants::over_c_sqr<double>);

      const auto& [
         mg1, mg2,
         gamma_cm, gamma1_cm, gamma2_cm,
         vcm_dot_v1, vcm_dot_v2, vcm,
         p1, p2, p1_star
      ] = getCOMData(particle1, particle2, g1.mass, g2.mass);

      const auto vcm_len2 = vcm.length_squared();
      const auto p1_star_len = p1_star.length();
      const auto gamma1_star = particle1.gamma * gamma1_cm;
      const auto gamma2_star = particle2.gamma * gamma2_cm;

      const auto gamma_coef1 = gamma_cm / (mg1 + mg2);
      const auto gamma_coef2 = 1.0 + constants::c_sqr<double> * (gamma1_star * g1.mass) * (gamma2_star * g2.mass) / math::SQR(p1_star_len);

      double coulomb_log_pairwise = coulombLog;
      if (coulombLog <= 0.0) {
         const auto l_deBroglie = h_half / p1_star_len;
         const auto b0 = gamma_coef1 * gamma_coef2 * coef2;
         const auto bmin2 = std::max(math::SQR(b0), math::SQR(l_deBroglie));
         coulomb_log_pairwise = std::max(2.0, 0.5 * std::log(1.0 + bmax2 / bmin2));
      }

      const auto s12_calc = rate_mult * max_weight * coef1 * scatter_coef * gamma_coef1 * math::SQR(gamma_coef2) * coulomb_log_pairwise * p1_star_len / (particle1.gamma * particle2.gamma * constants::c_sqr<double>);
      const auto s12_max = scatter_lowT * rate_mult * max_weight * dv_length;
      const auto s12 = std::min(s12_max, s12_calc);

      auto getScatteringAngles = [&]() -> std::array<double, 2> {
         double cos_theta;
         if (s12 >= 4.0) {
            cos_theta = 2.0 * rand1 - 1.0;
         }
         else {
            const auto alpha = 0.37 * s12 - 0.005 * math::SQR(s12) - 0.0064 * math::CUBE(s12);
            const auto sin2x2 = rand1 * alpha / std::sqrt(1.0 + rand1 * (math::SQR(alpha) - 1.0));
            cos_theta = 1.0 - (2.0 * sin2x2);
         }
         return {cos_theta, std::sqrt(1.0 - math::SQR(cos_theta))};
      };

      const auto [p1f_star, pcoef_lab] = calculate_P1final(p1_star, vcm, vcm_len2, gamma_cm, rand2, getScatteringAngles);
      applyScatter(particle1, particle2, vcm, p1f_star, pcoef_lab, gamma1_star, gamma2_star, gamma_cm, scatter_p1, scatter_p2);

      // const auto cos_theta = getCosTheta();
      // const auto sin_theta = std::sqrt(1.0 - math::SQR(cos_theta));
      // const auto phi = 2.0 * constants::pi<double> * rand2;
      // const auto dv0 = sin_theta * std::cos(phi);
      // const auto dv1 = sin_theta * std::sin(phi);
      //
      // const auto p1_star_perp = std::sqrt(math::SQR(p1_star[0]) + math::SQR(p1_star[1]));
      // const vec3 p1f_star{
      //    dv0 * (p1_star[0] * p1_star[2] / p1_star_perp) - dv1 * (p1_star[1] * p1_star_len / p1_star_perp) + cos_theta * p1_star[0],
      //    dv0 * (p1_star[1] * p1_star[2] / p1_star_perp) + dv1 * (p1_star[0] * p1_star_len / p1_star_perp) + cos_theta * p1_star[1],
      //    -dv0 * p1_star_perp + cos_theta * p1_star[2]
      // };
      //
      // const auto pcoef_lab = (gamma_cm - 1.0) * dot(vcm, p1f_star) / vcm_len2;
      // if (scatter_p1) {
      //    const auto p1f = p1f_star + (pcoef_lab + g1.mass * gamma1_star * gamma_cm) * vcm;
      //    const auto gamma1f = std::sqrt(1.0 + p1f.length_squared() * constants::over_c_sqr<double> / math::SQR(g1.mass));
      //    particle1.velocity = p1f / (g1.mass * gamma1f);
      //    particle1.gamma = gamma1f;
      // }
      //
      // if (scatter_p2) {
      //    const auto p2f = -p1f_star + (-pcoef_lab + g2.mass * gamma2_star * gamma_cm) * vcm;
      //    const auto gamma2f = std::sqrt(1.0 + p2f.length_squared() * constants::over_c_sqr<double> / math::SQR(g2.mass));
      //    particle2.velocity = p2f / (g2.mass * gamma2f);
      //    particle2.gamma = gamma2f;
      // }
   }

   void ionizationCollision(auto& particle1, auto& particle2,
      const auto weight1, const auto weight2, const auto max_weight,
      const auto dups, const auto scatter_coef,
      const auto rand0, const auto rand1, const auto rand2, const auto rand3, const auto rand4)
   {
      const auto dv_length = (particle1.velocity - particle2.velocity).length();
      const auto m_reduced = (g1.mass * g2.mass) / (g1.mass + g2.mass);

      // const auto mg1 = g1.mass * particle1.gamma;
      // const auto mg2 = g2.mass * particle2.gamma;
      //
      // const auto p1 = mg1 * particle1.velocity;
      // const auto p2 = mg2 * particle2.velocity;
      //
      // const auto vcm = (p1 + p2) / (mg1 + mg2);
      // const auto vcm_len2 = vcm.length_squared();
      // const auto gamma_cm = particles::calculateGammaV(vcm);
      //
      // const auto vcm_dot_v1 = dot(vcm, particle1.velocity);
      // const auto vcm_dot_v2 = dot(vcm, particle2.velocity);
      //
      // const auto gamma1_cm = gamma_cm * (particle1.gamma - vcm_dot_v1 * constants::over_c_sqr<double>);
      // const auto gamma2_cm = gamma_cm * (particle2.gamma - vcm_dot_v2 * constants::over_c_sqr<double>);

      const auto& [
         mg1, mg2,
         gamma_cm, gamma1_cm, gamma2_cm,
         vcm_dot_v1, vcm_dot_v2, vcm,
         p1, p2, p1_star
      ] = getCOMData(particle1, particle2, g1.mass, g2.mass);
      
      const auto vcm_len2 = vcm.length_squared();
      
      const vec3 v1_cm = particle1.velocity + vcm * ((gamma_cm - 1.0) / vcm_len2 * vcm_dot_v1 - gamma_cm * particle1.gamma);
      const vec3 v2_cm = particle2.velocity + vcm * ((gamma_cm - 1.0) / vcm_len2 * vcm_dot_v2 - gamma_cm * particle2.gamma);
      
      const auto dv2_cm = (v1_cm / gamma1_cm - v2_cm / gamma2_cm).length_squared();
      const auto vrel2_cm = dv2_cm / math::SQR(1.0 - dv2_cm * constants::over_c_sqr<double>);
      
      const auto gamma_rel = 1.0 / std::sqrt(1.0 - vrel2_cm * constants::over_c_sqr<double>);
      const auto energy_com_eV = (gamma_rel - 1.0) * m_reduced * constants::c_sqr<double> / constants::q_e<double>;

      if  (energy_com_eV < ionization_energy) { return; }

      const auto probability_coef = cross_section * max_weight * scatter_coef * dv_length * rate_mult;
      auto prod_mult = production_mult;
      auto probability_to_scatter = probability_coef * prod_mult;

      while (probability_to_scatter > probability_search_area and prod_mult > 1.0) {
         prod_mult /= 2.0;
         probability_to_scatter = probability_coef * prod_mult;
      }

      if (rand0 > probability_to_scatter) {
         return;
      }

      const auto Ee_scatter = rand1 * (energy_com_eV - ionization_energy);
      const auto pe_scatter = p1 * std::sqrt(Ee_scatter / energy_com_eV);
      const auto gamma_e_scatter = std::sqrt(1.0 + math::SQR(pe_scatter) * constants::over_c_sqr<double> / constants::m_e<double>);
      const auto v_e_scatter = pe_scatter / (gamma_e_scatter * g1.mass);

      const auto wi_eff = weight2 / rejection_mult;

      auto target_ionizes = weight1 >= wi_eff or (weight1 / wi_eff) >= rand2;
      auto electron_scatters = weight1 < wi_eff or (wi_eff / weight1) >= rand2;

      if (target_ionizes) {
         const auto product_weight = std::min(particle2.weight, wi_eff * dups / production_mult);
         
         const auto Ee_ejected = energy_com_eV - ionization_energy - Ee_scatter;
         const auto pe_ejected = p1 * std::sqrt(Ee_ejected / energy_com_eV);
         const auto gamma_e_ejected = std::sqrt(1.0 + math::SQR(pe_ejected) * constants::over_c_sqr<double> / constants::m_e<double>);
         const auto v_e_ejected = pe_ejected / (gamma_e_ejected * constants::m_e<double>);
         // todo: add new electron to relevant group here.
         product_electrons.add(
            Particle{
               particle2.location,
               particle2.location,
               v_e_ejected,
               gamma_e_ejected,
               product_weight,
               false
            }
         );
         
         // todo: add new ion to relevant group here.
         product_ions.add(
            Particle{
               particle2.location,
               particle2.location,
               particle2.velocity,
               particle2.gamma,
               product_weight,
               false
            }
         );

         if (product_weight == particle2.weight) {
            particle2.disabled = true;
         }
         else {
            particle2.weight -= product_weight;
         }

         // todo: accumulate ionization energy here
      } // end target_ionizes
      
      const auto gm1_gm2 = gamma_e_scatter * g1.mass + particle2.gamma * g2.mass;
      const auto vcm_scatter = (pe_scatter + p2) / gm1_gm2;
      const auto vcm_scatter_len2 = vcm_scatter.length_squared();

      const auto vcm_dot_ve_scatter = dot(vcm_scatter, v_e_scatter);
      const auto vcm_dot_v2_scatter = dot(vcm_scatter, particle2.velocity);
      const auto gamma_scatter_cm = particles::calculateGammaV(vcm_scatter);
      
      const auto pcoef = g1.mass * gamma_e_scatter * (vcm_dot_ve_scatter / vcm_scatter_len2 * (gamma_scatter_cm - 1.0) - gamma_scatter_cm);
      const auto p1_scatter_star = pe_scatter + pcoef * vcm_scatter;
      const auto p1_star_len = p1_scatter_star.length();

      const auto gamma1_star = gamma_e_scatter * gamma_scatter_cm * (1.0 - vcm_dot_ve_scatter * constants::over_c_sqr<double>);
      const auto gamma2_star = particle2.gamma * gamma_scatter_cm * (1.0 - vcm_dot_v2_scatter * constants::over_c_sqr<double>);

      auto getScatteringAngle = [&]() -> std::array<double, 2>
      {
         const auto& coefs = Park2022_ScatterCoefs[park_emission_target_type];
         const auto C = coefs[0] + coefs[1] * std::tanh(coefs[2] * (energy_com_eV - coefs[3]));
         const auto eta_f = coefs[4] / (energy_com_eV + coefs[5]) + coefs[6];
         const auto eta_b = coefs[7] / (energy_com_eV + coefs[8]) + coefs[9];

         const auto c1 = rand4 * eta_f + eta_b * (eta_f + 1.0 - rand4);
         const auto c2 = eta_f * eta_b + C * (eta_f + eta_b + 1.0) + rand4 * (eta_f - eta_b - 1.0);
         const auto c3 = eta_b + rand4 - C * (eta_f + eta_b + 1.0);
         const auto c4 = 4.0 * rand4 * eta_f * (eta_b + 1.0);

         const auto cos_theta = (c1 - std::sqrt(c2 * c2 + c3 * c4)) / c3;
         return {cos_theta, std::sqrt(1.0 - math::SQR(cos_theta))};
      };
      
      const auto [p1f_scatter_star, pcoef_scatter_lab] = calculate_P1final(p1_scatter_star, vcm_scatter, vcm_scatter_len2, gamma_scatter_cm, rand3, getScatteringAngle);
      applyScatter(particle1, particle2, vcm_scatter, p1f_scatter_star, pcoef_scatter_lab, gamma1_star, gamma2_star, gamma_scatter_cm, electron_scatters, target_ionizes);

      // const auto phi = 2.0 * constants::pi<double> * rand3;
      // // todo: simplify these
      // const auto theta = getScatteringAngle();
      // const auto dv0 = std::sin(theta) * std::cos(phi);
      // const auto dv1 = std::sin(theta) * std::sin(phi);
      // const auto cos_theta = std::cos(theta);
      //
      // const auto p1_star_perp = std::sqrt(math::SQR(p1_star[0]) + math::SQR(p1_star[1]));
      // const vec3 p1f_star{
      //    dv0 * (p1_star[0] * p1_star[2] / p1_star_perp) - dv1 * (p1_star[1] * p1_star_len / p1_star_perp) + cos_theta * p1_star[0],
      //    dv0 * (p1_star[1] * p1_star[2] / p1_star_perp) + dv1 * (p1_star[0] * p1_star_len / p1_star_perp) + cos_theta * p1_star[1],
      //    -dv0 * p1_star_perp + cos_theta * p1_star[2]
      // };
      //
      // const auto pcoef_lab = (gamma_cm - 1.0) * dot(vcm, p1f_star) / vcm_len2;

      // if (electron_scatters) {
      //    const auto p1f = p1f_star + (pcoef_lab + g1.mass * gamma1_star * gamma_cm) * vcm;
      //    const auto gamma1f = std::sqrt(1.0 + p1f.length_squared() * constants::over_c_sqr<double> / math::SQR(g1.mass));
      //    particle1.velocity = p1f / (g1.mass * gamma1f);
      //    particle1.gamma = gamma1f;
      // }
      //
      // if (target_ionizes) {
      //    const auto p2f = -p1f_star + (-pcoef_lab + g2.mass * gamma2_star * gamma_cm) * vcm;
      //    const auto gamma2f = std::sqrt(1.0 + p2f.length_squared() * constants::over_c_sqr<double> / math::SQR(g2.mass));
      //    particle2.velocity = p2f / (g2.mass * gamma2f);
      //    particle2.gamma = gamma2f;
      // }
   }

   void advance(const auto step) requires (coll_enabled) {
      static constexpr auto cell_vol_inv = 1.0 / (dx * dy * dz);
      static constexpr auto dt_vol = dt * cell_vol_inv;
      // static constexpr auto temp_coef  = 2.0 / (3.0 * constants::q_e<double>);
      static constexpr auto q_eps = constants::q_e<double> * constants::eps0<double>;
      static constexpr auto pi43 = 4.0 * constants::pi<double> / 3.0;
      static constexpr auto h_half = 0.5 * constants::h<double>;

      constexpr auto pi34_cuberoot = std::pow(pi43, 1.0 / 3.0);
      constexpr auto four_pi_eps_c2_inv = 1.0 / (4.0 * constants::pi<double> * constants::eps0<double> * constants::c_sqr<double>);

      if (step % step_interval != 0) { return; }

      const auto q1q2 = g1.charge * g2.charge;
      const auto m1m2 = g1.mass * g2.mass;
      const auto m1_over_m2 = g1.mass / g2.mass;
      const auto m1c2 = g1.charge * constants::c_sqr<double>;
      const auto m2c2 = g2.charge * constants::c_sqr<double>;
      const auto m_reduced = (g1.mass * g2.mass) / (g1.mass + g2.mass);

      const auto coef1 = four_pi_eps_c2_inv * math::SQR(q1q2) / (m1m2 * constants::eps0<double>);
      const auto coef2 = four_pi_eps_c2_inv * std::abs(q1q2);

      g1.update_cell_map();
      g2.update_cell_map();

      // pids1.clear();
      // pids2.clear();
      // nDups.clear();
      for (auto& [z_code, cell1] : g1.cell_map) {
         if (!g2.cell_map.contains(z_code)) {
            std::println("Cell {} not found in Group2.", z_code);
            continue;
         }

         const auto cell2 = g2.cell_map.at(z_code);
         const auto tid = 0; // omp_get_thread_num();

         const auto np1 = cell1.size();
         const auto np2 = cell2.size();

         const auto n_partners = self_scatter ? np1 - 1 + np1 % 2 : std::max(np1, np2);
         const auto scatter_coef = dt_vol * static_cast<double>(n_partners);

         std::uniform_real_distribution rng(0.0, 1.0);
         std::vector<std::size_t> pids1(n_partners);
         std::vector<std::size_t> pids2(n_partners);
         std::vector<double> nDups(n_partners);
         // pids1.resize(n_partners);
         // pids2.resize(n_partners);
         // nDups.resize(n_partners);

         auto calcDensityTemp = [&](const auto& c, const auto& mc2) -> std::tuple<double, double>
         {
            // Calculates cell density and temperature in eV
            auto ttl_weight = 0.0;
            auto KE = 0.0;
            for (auto& p : c) {
               ttl_weight += p.weight;
               KE += (p.gamma - 1.0) * p.weight * mc2;
            }
            return {ttl_weight * cell_vol_inv, (2.0 / 3.0) * KE / (ttl_weight * constants::q_e<double>)};
         };

         const auto& rhoT1 = calcDensityTemp(cell1, m1c2);
         const auto& [density1, temp1] = rhoT1;
         const auto& [density2, temp2] = self_scatter ? rhoT1 : calcDensityTemp(cell2, m2c2);

         // Cell Data
         const auto rmin2 = std::pow(pi43 * std::max(density1, density2), -2.0 / 3.0);
         const auto lD2_1 = q_eps * temp1 / (density1 * math::SQR(g1.charge));
         const auto lD2_2 = self_scatter ? 0.0 : q_eps * temp2 / (density2 * math::SQR(g2.charge));
         const auto bmax2 = rmin2 + lD2_1 + lD2_2;

         // todo: clean these lines up
         const auto d1 = m1_over_m2 * std::pow(density1, 2.0 / 3.0);
         const auto d2 = std::pow(density2, 2.0 / 3.0);
         const auto scatter_lowT = scatter_coef * pi34_cuberoot * (m1_over_m2 + 1.0) / std::max(d1, d2);

         // Create vector of random integers in range [0, Np1) and [0, Np2)
         for (std::size_t i = 0; i < n_partners; ++i) {
            pids1[i] = i % np1;
            pids2[i] = i % np2;
         }

         const auto grp_dups = static_cast<double>(np1 < np2 ? np2 / np1 : np1 / np2);
         const auto      rem = static_cast<int>(np1 < np2 ? np2 % np1 : np1 % np2);
         std::ranges::fill(nDups, grp_dups);
         if (rem != 0) {
            std::ranges::fill(nDups.end() - rem, nDups.end(), grp_dups);
         }

         std::ranges::shuffle(pids1, generator[tid]);
         std::ranges::shuffle(pids2, generator[tid]);

         #pragma omp parallel for num_threads(nThreads)
         for (std::size_t i = 0; i < n_partners; ++i) {
            auto& particle1 = cell1[pids1[i]];
            auto& particle2 = cell2[pids2[i]];

            const auto weight1 = particle1.weight / nDups[pids1[i]];
            const auto weight2 = particle2.weight / nDups[pids2[i]];
            const auto max_weight = std::max(weight1, weight2);

            coulombCollision(particle1, particle2, weight1, weight2, max_weight, rng(generator[tid]), rng(generator[tid]), rng(generator[tid]), coef1, coef2, bmax2, scatter_coef, scatter_lowT);


         } // end for(npairs)
      } // end for(z_code, cell1)

      // if (added_particles_to_group1) { g1.cell_map_updated = false; }
      // if (added_particles_to_group2) { g2.cell_map_updated = false; }

   } // end update()

   void advance(const auto) requires (!coll_enabled) {}

}; // end struct Collisions

} // end namespace tf::collisions

#endif //COLLISIONS_HPP
