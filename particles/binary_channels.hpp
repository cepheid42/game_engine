#ifndef GAME_ENGINE_BINARY_CHANNELS_HPP
#define GAME_ENGINE_BINARY_CHANNELS_HPP

#include "particles.hpp"
#include "constants.hpp"
#include "math_utils.hpp"
#include "rng_utils.h"
#include "vec3.hpp"

#include <algorithm>
#include <cmath>
#include <ranges>


namespace tf::collisions
{

struct ParticlePairData {
   particles::Particle& particle1;
   particles::Particle& particle2;
   double mass1;
   double mass2;
   double weight1;
   double weight2;
   double max_weight;
   double scatter_coef;
   std::array<double, 4> rand;
};
            
struct CoulombData {
   double coef1;
   double coef2;
   double bmax2;
   double scatter_lowt;
};

struct COMData {
   vec3<double> velocity;
   vec3<double> p1;
   vec3<double> p2;
   vec3<double> p1_star;
   double mg1;
   double mg2;
   double gamma;
   double gamma1;
   double gamma2;
   double vcm_dot_v1;
   double vcm_dot_v2;
   double vcm_len2;
};

auto calculate_P1final(const auto& p1_star, const auto& vcm, const auto vcm_len2, const auto gamma_cm, const auto rand0, const auto getScatteringAngle)
-> std::tuple<vec3<double>, double>
{
   const auto p1_star_len = p1_star.length();

   const auto [cos_theta, sin_theta] = getScatteringAngle();
   const auto phi = 2.0 * constants::pi * rand0;
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

auto updateParticles(auto& p1, auto& p2, const auto mass1, const auto mass2, const auto& vcm, const auto& p1f_star, const auto pcoef_lab, const auto g1_star, const auto g2_star, const auto g_cm, const auto cond1, const auto cond2) {
   if (cond1) {
      const auto p1f = p1f_star + (pcoef_lab + mass1 * g1_star * g_cm) * vcm;
      const auto gamma1f = particles::calculateGammaP(p1f, mass1);
      p1.velocity = p1f / (mass1 * gamma1f);
      p1.gamma = gamma1f;
   }

   if (cond2) {
      const auto p2f = -p1f_star + (-pcoef_lab + mass2 * g2_star * g_cm) * vcm;
      const auto gamma2f = particles::calculateGammaP(p2f, mass2);
      p2.velocity = p2f / (mass2 * gamma2f);
      p2.gamma = gamma2f;
   }
}

auto getCOMData(const auto& v1, const auto& v2, const auto mg1, const auto mg2)
-> COMData
{
   const auto p1 = mg1 * v1;
   const auto p2 = mg2 * v2;

   const auto vcm = (p1 + p2) / (mg1 + mg2);
   const auto vcm_len2 = vcm.length_squared();
   const auto gamma_cm = particles::calculateGammaV(vcm);

   const auto vcm_dot_v1 = dot(vcm, v1);
   const auto vcm_dot_v2 = dot(vcm, v2);

   const auto p1_star = p1 + vcm * mg1 * (vcm_dot_v1 / vcm_len2 * (gamma_cm - 1.0) - gamma_cm);

   const auto gamma1_cm = gamma_cm * (1.0 - vcm_dot_v1 * constants::over_c_sqr);
   const auto gamma2_cm = gamma_cm * (1.0 - vcm_dot_v2 * constants::over_c_sqr);

   return {
      .velocity = vcm,
      .p1 = p1,
      .p2 = p2,
      .p1_star = p1_star,
      .mg1 = mg1,
      .mg2 = mg2,
      .gamma = gamma_cm,
      .gamma1 = gamma1_cm,
      .gamma2 = gamma2_cm,
      .vcm_dot_v1 = vcm_dot_v1,
      .vcm_dot_v2 = vcm_dot_v2,
      .vcm_len2 = vcm.length_squared(),
   };
}

void coulombCollision(const bool has_coulomb, const auto& params, const auto& spec, const auto& coulomb)
{
   if (!has_coulomb) { return; } // no coulomb collisions

   auto& particle1 = params.particle1;
   auto& particle2 = params.particle2;
   const auto& cspec = spec.coulomb;

   const auto scatter_p1 = params.rand[0] * params.weight1 < params.weight2;
   const auto scatter_p2 = params.rand[0] * params.weight2 < params.weight1;
   const auto dv_length = (particle1.velocity - particle2.velocity).length();

   if (dv_length == 0.0 or !(scatter_p1 or scatter_p2)) {
      return;
   } // parallel particles or self particles can cause divide by zero errors

   const auto com = getCOMData(particle1.velocity, particle2.velocity, particle1.gamma * params.mass1, particle2.gamma * params.mass2);

   const auto p1_star_len = com.p1_star.length();
   const auto gamma1_star = particle1.gamma * com.gamma1;
   const auto gamma2_star = particle2.gamma * com.gamma2;

   const auto gamma_coef1 = com.gamma / (com.mg1 + com.mg2);
   const auto gamma_coef2 = 1.0 + constants::c_sqr * (gamma1_star * params.mass1) * (gamma2_star * params.mass2) / math::SQR(p1_star_len);

   double coulomb_log_pairwise = cspec.coulomb_log;
   if (coulomb.log <= 0.0) {
      const auto l_deBroglie = 0.5 * constants::h / p1_star_len;
      const auto b0 = gamma_coef1 * gamma_coef2 * coulomb.coef2;
      const auto bmin2 = std::max(math::SQR(b0), math::SQR(l_deBroglie));
      coulomb_log_pairwise = std::max(2.0, 0.5 * std::log(1.0 + coulomb.bmax2 / bmin2));
   }

   const auto s12_calc = cspec.rate_mult * params.max_weight * coulomb.coef1 * params.scatter_coef * gamma_coef1 * math::SQR(gamma_coef2) * coulomb_log_pairwise * p1_star_len / (particle1.gamma * particle2.gamma * constants::c_sqr);
   const auto s12_max = coulomb.scatter_lowt * cspec.rate_multplier * params.max_weight * dv_length;
   const auto s12 = std::min(s12_max, s12_calc);

   auto getScatteringAngles = [&]() -> std::array<double, 2> {
      double cos_theta;
      if (s12 >= 4.0) {
         cos_theta = 2.0 * params.rand[1] - 1.0;
      }
      else {
         const auto alpha = 0.37 * s12 - 0.005 * math::SQR(s12) - 0.0064 * math::CUBE(s12);
         const auto sin2x2 = params.rand[1] * alpha / std::sqrt(1.0 + params.rand[1] * (math::SQR(alpha) - 1.0));
         cos_theta = 1.0 - (2.0 * sin2x2);
      }
      return {cos_theta, std::sqrt(1.0 - math::SQR(cos_theta))};
   };

   const auto& [p1f_star, pcoef_lab] = calculate_P1final(com.p1_star, com.velocity, com.vcm_len2, com.gamma, params.rand[2], getScatteringAngles);
   updateParticles(
      particle1, particle2,
      com.velocity, p1f_star, pcoef_lab,
      gamma1_star, gamma2_star, com.gamma,
      scatter_p1, scatter_p2
   );
}

void ionizationCollision(const bool has_ionization, const auto& params, const auto& spec, auto& product1, auto& product2, const auto ndups)
{
   if (!has_ionization) { return; } // no ionization collisions

   auto& particle1 = params.particle1;
   auto& particle2 = params.particle2;
   const auto& ionization = spec.ionization;
   
   const auto wi_eff = params.weight2 / ionization.rejection_mult;
   const auto dv_length = (particle1.velocity - particle2.velocity).length();
   const auto m_reduced = (params.mass1 * params.mass2) / (params.mass1 + params.mass2);

   auto target_ionizes = params.weight1 >= wi_eff or (params.weight1 / wi_eff >= params.rand[2]);
   auto electron_scatters = params.weight1 < wi_eff or (wi_eff / params.weight1 >= params.rand[2]);


   if (dv_length == 0.0 or !(target_ionizes or electron_scatters)) {
      return;
   } // parallel particles or self particles can cause divide by zero errors

   const auto& com = getCOMData(particle1.velocity, particle2.velocity, particle1.gamma * params.mass1, particle2.gamma * params.mass2);

   const vec3 v1_cm = particle1.velocity + com.velocity * ((com.gamma - 1.0) / com.vcm_len2 * com.vcm_dot_v1 - com.gamma * particle1.gamma);
   const vec3 v2_cm = particle2.velocity + com.velocity * ((com.gamma - 1.0) / com.vcm_len2 * com.vcm_dot_v2 - com.gamma * particle2.gamma);

   const auto dv2_cm = (v1_cm / com.gamma1 - v2_cm / com.gamma2).length_squared();
   const auto vrel2_cm = dv2_cm * particles::calculateGammaV(dv2_cm);

   const auto gamma_rel = calculateGamma(vrel2_cm);
   const auto energy_com_eV = (gamma_rel - 1.0) * m_reduced * constants::c_sqr / constants::q_e;

   double cross_section{ionization.constant_cross_section};
   if (cross_section == 0.0) {
      // todo: fill this in with table lookup
   }

   const auto probability_coef = cross_section * params.max_weight * params.scatter_coef * dv_length * ionization.rate_mult;
   auto prod_mult = ionization.production_mult;
   auto scatter_probability = probability_coef * prod_mult;

   while (scatter_probability > spec.search_area_probability and prod_mult > 1.0) {
      prod_mult /= 2.0;
      scatter_probability = probability_coef * prod_mult;
   }

   // Skip if COM energy is below threshold or scatter probability is below random number
   if (energy_com_eV < ionization.ionization_energy or params.rand[0] > scatter_probability) { return; }

   const auto Ee_scatter = params.rand[1] * (energy_com_eV - ionization.ionization_energy);
   const auto pe_scatter = com.p1 * std::sqrt(Ee_scatter / energy_com_eV);
   const auto gamma_e_scatter = calculateGammaP(pe_scatter, constants::m_e);
   const auto ve_scatter = pe_scatter / (gamma_e_scatter * params.mass1);

   if (target_ionizes) {
      const auto product_weight = std::min(particle2.weight, wi_eff * ndups / ionization.production_mult);
      const auto Ee_ejected = energy_com_eV - ionization.ionization_energy - Ee_scatter;
      const auto pe_ejected = com.p1 * std::sqrt(Ee_ejected / energy_com_eV);
      const auto gamma_e_ejected = particles::calculateGammaP(pe_ejected, constants::m_e);
      const auto v_e_ejected = pe_ejected / (gamma_e_ejected * constants::m_e);

      // todo: is this enough to keep OpenMP threads from stepping on each other?
      #pragma omp critical
      {
         product1.emplace_back(
            v_e_ejected,
            gamma_e_ejected,
            particle2.location,
            particle2.location,
            product_weight
         );

         product2.emplace_back(
            particle2.velocity,
            particle2.gamma,
            particle2.location,
            particle2.location,
            product_weight
         );
      }
      if (product_weight == particle2.weight) {
         particle2.weight = -1.0;
      }
      else {
         particle2.weight -= product_weight;
      }
   } // end target_ionizes

   // particles::Particle e_scattered{{}, {}, ve_scatter, 0.0, gamma_e_scatter, false};

   const auto& com_scatter = getCOMData(ve_scatter, particle2.velocity, gamma_e_scatter * params.mass1, particle2.gamma * params.mass2);

   const auto gamma1_star = gamma_e_scatter * com_scatter.gamma * (1.0 - com_scatter.vcm_dot_ve * constants::over_c_sqr);
   const auto gamma2_star = particle2.gamma * com_scatter.gamma * (1.0 - com_scatter.vcm_dot_v2 * constants::over_c_sqr);

   auto getScatteringAngle = [&]() -> std::array<double, 2> {
      const auto cos_theta = 2.0 * params.rand[1] - 1.0;
      return {cos_theta, std::sqrt(1.0 - math::SQR(cos_theta))};
   };

   const auto& [p1f_star_scatter, pcoef_lab_scatter] = calculate_P1final(com_scatter.p1_star, com_scatter.velocity, com_scatter.vcm_len2, com_scatter.gamma_cm, params.rand[3], getScatteringAngle);
   updateParticles(
      particle1, particle2,
      com_scatter.velocity, p1f_star_scatter, pcoef_lab_scatter,
      gamma1_star, gamma2_star, com_scatter.gamma,
      electron_scatters, target_ionizes
   );
}
} // end namespace tf::collisions

#endif //GAME_ENGINE_BINARY_CHANNELS_HPP