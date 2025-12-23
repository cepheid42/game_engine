#ifndef GAME_ENGINE_BINARY_CHANNELS_HPP
#define GAME_ENGINE_BINARY_CHANNELS_HPP

#include "constants.hpp"
#include "math_utils.hpp"
#include "particles.hpp"
#include "rng_utils.h"
#include "vec3.hpp"

#include <algorithm>
#include <cmath>
#include <ranges>
#include <cassert>

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
   std::array<double, 5> rand;
};
            
struct CoulombData {
   double coef1;
   double coef2;
   double bmax2;
   double scatter_lowt;
};

auto computecCOMEnergy(auto m1, auto m2, auto particle1, auto particle2) {
   const auto mg1 = m1 * particle1.gamma;
   const auto mg2 = m2 * particle2.gamma;

   const auto mu = (m1 * m2) / (m1 + m2);

   const auto p1 = mg1 * particle1.velocity;
   const auto p2 = mg2 * particle2.velocity;

   // std::println("momentum: {} ||| {}", p1, p2);

   const auto vcm = (p1 + p2) / (mg1 + mg2);
   const auto vcm_sqr = vcm.length_squared();

   // std::println("vcm_sqr: {}", vcm_sqr);

   const auto gamma_cm = 1.0 / std::sqrt(1.0 - vcm_sqr * constants::over_c_sqr);

   // std::println("gamma_cm: {}", gamma_cm);

   const auto u1 = p1 / m1;
   const auto u2 = p2 / m2;

   const auto vcm_dot_u1 = dot(vcm, u1);
   const auto vcm_dot_u2 = dot(vcm, u2);

   // std::println("vcm_dot_u1: {}", vcm_dot_u1);
   // std::println("vcm_dot_u2: {}", vcm_dot_u2);

   const auto gamma1_cm = gamma_cm * (particle1.gamma - vcm_dot_u1 * constants::over_c_sqr);
   const auto gamma2_cm = gamma_cm * (particle2.gamma - vcm_dot_u2 * constants::over_c_sqr);

   // std::println("gamma12_cm: {} ||| {}", gamma1_cm, gamma2_cm);

   const auto u1_cm = u1 + vcm * ((gamma_cm - 1.0) * vcm_dot_u1 / vcm_sqr - gamma_cm * particle1.gamma);
   const auto u2_cm = u2 + vcm * ((gamma_cm - 1.0) * vcm_dot_u2 / vcm_sqr - gamma_cm * particle2.gamma);

   // std::println("u1_cm: {}", u1_cm);
   // std::println("u2_cm: {}", u2_cm);

   const auto dv2_cm = (u1_cm / gamma1_cm - u2_cm / gamma2_cm).length_squared();
   const auto vrel2_cm = dv2_cm / math::SQR(1.0 - dv2_cm * constants::over_c_sqr);

   const auto gamma_rel = 1.0 / std::sqrt(1.0 - vrel2_cm * constants::over_c_sqr);

   // std::println("{}", (gamma_rel - 1.0) * mu * constants::c_sqr / constants::q_e);

   return (gamma_rel - 1.0) * mu * constants::c_sqr / constants::q_e;
}

void ionizationCollision(
   const bool has_ionization,
   const auto& params,
   const auto& spec,
   auto& product1, auto& product2,
   const auto ndups, const auto& cs_table)
{
   if (!has_ionization) { return; } // no ionization collisions

   auto& particle1 = params.particle1;
   auto& particle2 = params.particle2;
   const auto& ionization = spec.ionization;

   const auto energy_com_eV = computecCOMEnergy(params.mass1, params.mass2, particle1, particle2);

   if (energy_com_eV <= ionization.ionization_energy) { return; }

   auto cross_section{ionization.constant_cross_section};
   if (cross_section == 0.0) {
      cross_section = cs_table.lerp(energy_com_eV);
   }

   const auto dv_length = (particle1.velocity - particle2.velocity).length();

   auto prod_mult = ionization.production_multiplier;
   const auto probability_coef = cross_section * params.max_weight * params.scatter_coef * dv_length * ionization.rate_multiplier;
   auto probability_to_scatter = probability_coef * prod_mult;
   while (probability_to_scatter > spec.probability_search_area and prod_mult > 1.0) {
      prod_mult /= 2.0;
      probability_to_scatter = probability_coef * prod_mult;
   }

   if (probability_to_scatter >= params.rand[0]) { return; }

   const auto p1 = params.mass1 * particle1.gamma * particle1.velocity;
   const auto p2 = params.mass2 * particle2.gamma * particle2.velocity;

   const auto E_scattered_electron = params.rand[1] * (energy_com_eV - ionization.ionization_energy);
   const auto p_e_scattered = p1 * std::sqrt(E_scattered_electron / energy_com_eV);
   const auto gamma_e_scattered = particles::calculateGammaP(p_e_scattered, constants::m_e);
   const auto v_e_scattered = p_e_scattered / (gamma_e_scattered * params.mass1);

   const auto E_ejected_electron = energy_com_eV - (ionization.ionization_energy + E_scattered_electron);
   const auto p_e_ejected = p1 * std::sqrt(E_ejected_electron / energy_com_eV);
   const auto gamma_e_ejected = particles::calculateGammaP(p_e_ejected, constants::m_e);
   const auto v_e_ejected = p_e_ejected / (gamma_e_ejected * params.mass1);

   const auto wi_effective = params.weight2 / ionization.rejection_multiplier;

   auto target_ionizes = params.weight1 >= wi_effective or (params.weight1 / wi_effective >= params.rand[2]);
   auto electron_scatters = params.weight1 < wi_effective or (wi_effective / params.weight1 >= params.rand[2]);

   if (target_ionizes) {
      const auto product_weight = std::min(particle2.weight, wi_effective * ndups / ionization.production_multiplier);

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
   }

   const auto gm1_gm2 = gamma_e_scattered * params.mass1 + particle2.gamma * params.mass2;
   const auto vcm_scattered = (p_e_scattered + p2) / gm1_gm2;
   const auto vcm_scattered_sqr = vcm_scattered.length_squared();

   const auto gamma_cm_scattered = particles::calculateGammaV(vcm_scattered);
   const auto pcoef = params.mass1 * gamma_e_scattered * ((gamma_cm_scattered - 1.0) * dot(vcm_scattered, v_e_scattered) / vcm_scattered_sqr - gamma_cm_scattered);
   const auto p1_star = p_e_scattered + pcoef * vcm_scattered;

   const auto p1_star_sqr = p1_star.length_squared();
   const auto p1_star_len = std::sqrt(p1_star_sqr);

   const auto p1_star_perp_sqr = math::SQR(p1_star[0]) + math::SQR(p1_star[1]);
   const auto p1_star_perp = std::sqrt(p1_star_perp_sqr);

   if (p1_star_sqr <= 0.0 or p1_star_perp_sqr <= 0.0) { return; }

   const auto gamma1_star = (1.0 - dot(vcm_scattered, v_e_scattered) * constants::over_c_sqr) * gamma_e_scattered * gamma_cm_scattered;
   const auto gamma2_star = (1.0 - dot(vcm_scattered, particle2.velocity) * constants::over_c_sqr) * particle2.gamma * gamma_cm_scattered;

   const auto phi = 2.0 * constants::pi * params.rand[3];

   auto getScatteringAngle = [&]() -> std::array<double, 2> {
      const auto cos_theta = 2.0 * params.rand[4] - 1.0;
      return {cos_theta, std::sqrt(1.0 - math::SQR(cos_theta))};
   };

   const auto [cos_theta, sin_theta] = getScatteringAngle();

   const auto dv0 = sin_theta * std::cos(phi);
   const auto dv1 = sin_theta * std::sin(phi);
   const auto dv2 = cos_theta;

   const auto p1f_star = vec3{
      dv0 * (p1_star[0] * p1_star[2] / p1_star_perp) - dv1 * (p1_star[1] * p1_star_len / p1_star_perp) + dv2 * p1_star[0],
      dv0 * (p1_star[1] * p1_star[2] / p1_star_perp) + dv1 * (p1_star[0] * p1_star_len / p1_star_perp) + dv2 * p1_star[1],
      -dv0 * p1_star_perp + dv2 * p1_star[2]
    };

   const auto pcoef_lab = (gamma_cm_scattered - 1.0) * dot(vcm_scattered, p1f_star) / vcm_scattered_sqr;
   const auto p1f = p1f_star + vcm_scattered * (pcoef_lab + params.mass1 * gamma1_star * gamma_cm_scattered);
   const auto p2f = -p1f_star + vcm_scattered * (-pcoef_lab + params.mass2 * gamma2_star * gamma_cm_scattered);

   const auto gamma1f = particles::calculateGammaP(p1f, params.mass1);
   const auto gamma2f = particles::calculateGammaP(p2f, params.mass2);

   if (electron_scatters) {
      particle1.velocity = p1f / (params.mass1 * gamma1f);
      particle1.gamma = gamma1f;
   }

   if (target_ionizes) {
      particle2.velocity = p2f / (params.mass2 * gamma2f);
      particle2.gamma = gamma2f;
   }
}

// struct COMData {
//    vec3<double> velocity;
//    vec3<double> p1;
//    vec3<double> p2;
//    vec3<double> p1_star;
//    double mg1;
//    double mg2;
//    double gamma;
//    double gamma1;
//    double gamma2;
//    double vcm_dot_v1;
//    double vcm_dot_v2;
//    double vcm_len2;
// };
//
// auto calculate_P1final(const auto& p1_star, const auto& vcm, const auto vcm_len2, const auto gamma_cm, const auto rand, const auto getScatteringAngle)
// -> std::array<vec3<double>, 2>
// {
//    const auto p1_star_len = p1_star.length();
//
//    const auto [cos_theta, sin_theta] = getScatteringAngle();
//    const auto phi = 2.0 * constants::pi * rand;
//    const auto dv0 = sin_theta * std::cos(phi);
//    const auto dv1 = sin_theta * std::sin(phi);
//
//    // Perez (2012) eq 12
//    const auto p1_star_perp = std::sqrt(math::SQR(p1_star[0]) + math::SQR(p1_star[1]));
//    const vec3 p1f_star {
//       dv0 * (p1_star[0] * p1_star[2] / p1_star_perp) - dv1 * (p1_star[1] * p1_star_len / p1_star_perp) + cos_theta * p1_star[0],
//       dv0 * (p1_star[1] * p1_star[2] / p1_star_perp) + dv1 * (p1_star[0] * p1_star_len / p1_star_perp) + cos_theta * p1_star[1],
//       -dv0 * p1_star_perp + cos_theta * p1_star[2]
//    };
//
//    const auto pcoef_lab = (gamma_cm - 1.0) * dot(vcm, p1f_star) / vcm_len2;
//
//    return {p1f_star, pcoef_lab};
// }
//
// auto updateParticles(auto& p1, auto& p2, const auto m1, const auto m2, const auto& vcm, const auto& p1f_star, const auto pcoef_lab, const auto g1_star, const auto g2_star, const auto g_cm, const auto cond1, const auto cond2) {
//    // Perez (2012) eq 13
//    if (cond1) {
//       const auto p1f = p1f_star + vcm * (pcoef_lab + m1 * g1_star * g_cm);
//       const auto gamma1f = particles::calculateGammaP(p1f, m1);
//       p1.velocity = p1f / (m1 * gamma1f);
//       p1.gamma = gamma1f;
//    }
//
//    if (cond2) {
//       const auto p2f = -p1f_star + vcm * (-pcoef_lab + m2 * g2_star * g_cm);
//       const auto gamma2f = particles::calculateGammaP(p2f, m2);
//       p2.velocity = p2f / (m2 * gamma2f);
//       p2.gamma = gamma2f;
//    }
// }
//
// auto getCOMData(const auto& v1, const auto& v2, const auto g1, const auto g2, const auto m1, const auto m2)
// -> COMData
// {
//    const auto mg1 = m1 * g1;
//    const auto mg2 = m2 * g2;
//
//    const auto p1 = mg1 * v1;
//    const auto p2 = mg2 * v2;
//
//    const auto vcm = (p1 + p2) / (mg1 + mg2);
//    const auto vcm_len2 = vcm.length_squared();
//    const auto gamma_cm = particles::calculateGammaV(vcm);
//
//    const auto vcm_dot_v1 = dot(vcm, v1);
//    const auto vcm_dot_v2 = dot(vcm, v2);
//
//    // Perez (2012) eq 2
//    const auto p1_star = p1 + vcm * mg1 * ((gamma_cm - 1.0) * vcm_dot_v1 / vcm_len2 - gamma_cm);
//
//    const auto gamma1_cm = g1 * gamma_cm * (1.0 - vcm_dot_v1 * constants::over_c_sqr);
//    const auto gamma2_cm = g2 * gamma_cm * (1.0 - vcm_dot_v2 * constants::over_c_sqr);
//
//    return {
//       .velocity = vcm,
//       .p1 = p1,
//       .p2 = p2,
//       .p1_star = p1_star,
//       .mg1 = mg1,
//       .mg2 = mg2,
//       .gamma = gamma_cm,
//       .gamma1 = gamma1_cm,
//       .gamma2 = gamma2_cm,
//       .vcm_dot_v1 = vcm_dot_v1,
//       .vcm_dot_v2 = vcm_dot_v2,
//       .vcm_len2 = vcm.length_squared(),
//    };
// }
//
// auto calculateEnergyCOM(const auto& p1, const auto& p2, const auto& com, const auto m_reduced) {
//    const vec3 u1_cm = p1.gamma * p1.velocity + com.velocity * ((com.gamma - 1.0) * com.vcm_dot_v1 / com.vcm_len2 - com.gamma * p1.gamma);
//    const vec3 u2_cm = p2.gamma * p2.velocity + com.velocity * ((com.gamma - 1.0) * com.vcm_dot_v2 / com.vcm_len2 - com.gamma * p2.gamma);
//    const auto dv2_cm = (u1_cm / com.gamma1 - u2_cm / com.gamma2).length_squared();
//    const auto vrel2_cm = dv2_cm / math::SQR(1.0 - dv2_cm * constants::over_c_sqr);
//    const auto gamma_rel = 1.0 / std::sqrt(1.0 - vrel2_cm * constants::over_c_sqr);
//
//    return (gamma_rel - 1.0) * m_reduced * constants::c_sqr / constants::q_e;
// }
//
// auto computeCOMEnergy(auto m1, auto m2, auto particle1, auto particle2) {
//    const auto mg1 = m1 * particle1.gamma;
//    const auto mg2 = m2 * particle2.gamma;
//
//    const auto mu = (m1 * m2) / (m1 + m2);
//
//    const auto p1 = mg1 * particle1.velocity;
//    const auto p2 = mg2 * particle2.velocity;
//
//    const auto vcm = (p1 + p2) / (mg1 + mg2);
//    const auto vcm_sqr = vcm.length_squared();
//
//    const auto gamma_cm = 1.0 / std::sqrt(1.0 - vcm_sqr * constants::over_c_sqr);
//
//    const auto u1 = p1 / m1;
//    const auto u2 = p2 / m2;
//
//    const auto vcm_dot_u1 = dot(vcm, u1);
//    const auto vcm_dot_u2 = dot(vcm, u2);
//
//    const auto gamma1_cm = gamma_cm * (particle1.gamma - vcm_dot_u1 * constants::over_c_sqr);
//    const auto gamma2_cm = gamma_cm * (particle2.gamma - vcm_dot_u2 * constants::over_c_sqr);
//
//    const auto u1_cm = u1 + vcm * ((gamma_cm - 1.0) * vcm_dot_u1 / vcm_sqr - gamma_cm * particle1.gamma);
//    const auto u2_cm = u2 + vcm * ((gamma_cm - 1.0) * vcm_dot_u2 / vcm_sqr - gamma_cm * particle2.gamma);
//
//    const auto dv2_cm = (u1_cm / gamma1_cm - u2_cm / gamma2_cm).length_squared();
//    const auto vrel2_cm = dv2_cm / math::SQR(1.0 - dv2_cm * constants::over_c_sqr);
//
//    const auto gamma_rel = 1.0 / std::sqrt(1.0 - vrel2_cm * constants::over_c_sqr);
//
//    return (gamma_rel - 1.0) * mu * constants::c_sqr / constants::q_e;
// }
//
// void ionizationCollision(const bool has_ionization,
//    const auto& params,
//    const auto& spec,
//    auto& product1, auto& product2,
//    const auto ndups, const auto& cs_table)
// {
//    if (!has_ionization) { return; } // no ionization collisions
//
//    auto& particle1 = params.particle1;
//    auto& particle2 = params.particle2;
//    const auto& ionization = spec.ionization;
//
//    // const auto wi_eff = params.weight2 / ionization.rejection_multiplier;
//    // const auto dv_length = (particle1.velocity - particle2.velocity).length();
//    const auto m_reduced = (params.mass1 * params.mass2) / (params.mass1 + params.mass2);
//
//    // auto target_ionizes = params.weight1 >= wi_eff or (params.weight1 / wi_eff >= params.rand[2]);
//    // auto electron_scatters = params.weight1 < wi_eff or (wi_eff / params.weight1 >= params.rand[2]);
//    //
//    // if (dv_length == 0.0 or !(target_ionizes or electron_scatters)) {
//    //    return;
//    // } // parallel particles or self particles can cause divide by zero errors
//
//    // const auto com = getCOMData(particle1.velocity, particle2.velocity, particle1.gamma, particle2.gamma, params.mass1, params.mass2);
//    // const auto energy_com_eV = calculateEnergyCOM(particle1, particle2, com, m_reduced);
//
//    const auto mg1_a = params.mass1 * particle1.gamma;
//    const auto mg2_a = params.mass2 * particle2.gamma;
//
//    const auto p1_a = mg1_a * particle1.velocity;
//    const auto p2_a = mg2_a * particle2.velocity;
//
//    const auto vcm_a = (p1_a + p2_a) / (mg1_a + mg2_a);
//    const auto vcm_len2_a = vcm_a.length_squared();
//    const auto gamma_cm_a = particles::calculateGammaV(vcm_a);
//
//    const auto vcm_dot_v1_a = dot(vcm_a, particle1.velocity);
//    const auto vcm_dot_v2_a = dot(vcm_a, particle2.velocity);
//
//    const auto gamma1_cm_a = particle1.gamma * gamma_cm_a * (1.0 - vcm_dot_v1_a * constants::over_c_sqr);
//    const auto gamma2_cm_a = particle2.gamma * gamma_cm_a * (1.0 - vcm_dot_v2_a * constants::over_c_sqr);
//
//    const vec3 u1_cm_a = particle1.gamma * particle1.velocity + vcm_a * ((gamma_cm_a - 1.0) * vcm_dot_v1_a / vcm_len2_a - gamma_cm_a * particle1.gamma);
//    const vec3 u2_cm_a = particle2.gamma * particle2.velocity + vcm_a * ((gamma_cm_a - 1.0) * vcm_dot_v2_a / vcm_len2_a - gamma_cm_a * particle2.gamma);
//
//    const auto dv2_cm_a = (u1_cm_a / gamma1_cm_a - u2_cm_a / gamma2_cm_a).length_squared();
//    const auto vrel2_cm_a = dv2_cm_a / math::SQR(1.0 - dv2_cm_a * constants::over_c_sqr);
//    const auto gamma_rel_a = 1.0 / std::sqrt(1.0 - vrel2_cm_a * constants::over_c_sqr);
//
//    const auto energy_a = (gamma_rel_a - 1.0) * m_reduced * constants::c_sqr / constants::q_e;
//
//    // -----------------------------------------------------------------------------------------------
//    // const auto energy = computeCOMEnergy(params.mass1, params.mass2, particle1, particle2);
//
//    const auto mg1_b = params.mass1 * particle1.gamma;
//    const auto mg2_b = params.mass2 * particle2.gamma;
//
//    std::println("mg: {}, {}", mg1_a, mg1_b);
//
//    const auto p1_b = mg1_b * particle1.velocity;
//    const auto p2_b = mg2_b * particle2.velocity;
//
//    std::println("momentum: {}, {} |||| {}, {}", p1_a, p1_b, p2_a, p2_b);
//
//    const auto vcm_b = (p1_b + p2_b) / (mg1_b + mg2_b);
//    const auto vcm_sqr_b = vcm_b.length_squared();
//
//    std::println("vcm_sqr: {}, {}", vcm_len2_a, vcm_sqr_b);
//
//    const auto gamma_cm_b = 1.0 / std::sqrt(1.0 - vcm_sqr_b * constants::over_c_sqr);
//
//    std::println("gamma_cm: {}, {}", gamma_cm_a, gamma_cm_b);
//
//    const auto u1_b = p1_b / params.mass1;
//    const auto u2_b = p2_b / params.mass2;
//
//    const auto vcm_dot_u1_b = dot(vcm_b, u1_b);
//    const auto vcm_dot_u2_b = dot(vcm_b, u2_b);
//
//    std::println("vcm_dot_u1: {}, {}", vcm_dot_v1_a, vcm_dot_u1_b);
//    std::println("vcm_dot_u2: {}, {}", vcm_dot_v2_a, vcm_dot_u2_b);
//
//    const auto gamma1_cm_b = gamma_cm_b * (particle1.gamma - vcm_dot_u1_b * constants::over_c_sqr);
//    const auto gamma2_cm_b = gamma_cm_b * (particle2.gamma - vcm_dot_u2_b * constants::over_c_sqr);
//
//    std::println("gamma12_cm: {}, {} ||| {}, {}", gamma1_cm_a, gamma1_cm_b, gamma2_cm_a, gamma2_cm_b);
//
//    const auto u1_cm_b = u1_b + vcm_b * ((gamma_cm_b - 1.0) * vcm_dot_u1_b / vcm_sqr_b - gamma_cm_b * particle1.gamma);
//    const auto u2_cm_b = u2_b + vcm_b * ((gamma_cm_b - 1.0) * vcm_dot_u2_b / vcm_sqr_b - gamma_cm_b * particle2.gamma);
//
//    std::println("u1_cm: {}, {}", u1_cm_a, u1_cm_b);
//    std::println("u2_cm: {}, {}", u2_cm_a, u2_cm_b);
//
//    const auto dv2_cm_b = (u1_cm_b / gamma1_cm_b - u2_cm_b / gamma2_cm_b).length_squared();
//    const auto vrel2_cm_b = dv2_cm_b / math::SQR(1.0 - dv2_cm_b * constants::over_c_sqr);
//
//    const auto gamma_rel_b = 1.0 / std::sqrt(1.0 - vrel2_cm_b * constants::over_c_sqr);
//
//    const auto energy_b = (gamma_rel_b - 1.0) * m_reduced * constants::c_sqr / constants::q_e;
//
//    std::println("{}, {}", energy_a, energy_b);
//
//    // auto cross_section{ionization.constant_cross_section};
//    // if (cross_section == 0.0) {
//    //    cross_section = cs_table.lerp(energy_com_eV);
//    // }
//    //
//    // const auto probability_coef = cross_section * params.max_weight * params.scatter_coef * dv_length * ionization.rate_multiplier;
//    // auto prod_mult = ionization.production_multiplier;
//    // auto scatter_probability = probability_coef * prod_mult;
//    //
//    // while (scatter_probability > spec.probability_search_area and prod_mult > 1.0) {
//    //    prod_mult /= 2.0;
//    //    scatter_probability = probability_coef * prod_mult;
//    // }
//    //
//    // // Skip if COM energy is below threshold or scatter probability is below random number
//    // if (energy_com_eV <= ionization.ionization_energy or params.rand[0] > scatter_probability) { return; }
//    //
//    // const auto Ee_scatter = params.rand[1] * (energy_com_eV - ionization.ionization_energy);
//    // const auto pe_scatter = particle1.velocity * particle1.gamma * params.mass1 * std::sqrt(Ee_scatter / energy_com_eV);
//    // const auto gamma_e_scatter = particles::calculateGammaP(pe_scatter, constants::m_e);
//    // const auto ve_scatter = pe_scatter / (gamma_e_scatter * params.mass1);
//    //
//    // if (target_ionizes) {
//    //    const auto product_weight = std::min(particle2.weight, wi_eff * ndups / ionization.production_multiplier);
//    //    const auto Ee_ejected = energy_com_eV - (ionization.ionization_energy + Ee_scatter);
//    //    const auto pe_ejected = particle1.velocity * particle1.gamma * params.mass1 * std::sqrt(Ee_ejected / energy_com_eV);
//    //    const auto gamma_e_ejected = particles::calculateGammaP(pe_ejected, constants::m_e);
//    //    const auto v_e_ejected = pe_ejected / (gamma_e_ejected * constants::m_e);
//    //
//    //    // todo: is this enough to keep OpenMP threads from stepping on each other?
//    //    #pragma omp critical
//    //    {
//    //       product1.emplace_back(
//    //          v_e_ejected,
//    //          gamma_e_ejected,
//    //          particle2.location,
//    //          particle2.location,
//    //          product_weight
//    //       );
//    //
//    //       product2.emplace_back(
//    //          particle2.velocity,
//    //          particle2.gamma,
//    //          particle2.location,
//    //          particle2.location,
//    //          product_weight
//    //       );
//    //    }
//    //    if (product_weight == particle2.weight) {
//    //       particle2.weight = -1.0;
//    //    }
//    //    else {
//    //       particle2.weight -= product_weight;
//    //    }
//    // } // end target_ionizes
//    //
//    // const auto com_scatter = getCOMData(ve_scatter, particle2.velocity, gamma_e_scatter, particle2.gamma, params.mass1, params.mass2);
//    //
//    // const auto gamma1_star = gamma_e_scatter * com_scatter.gamma * (1.0 - com_scatter.vcm_dot_v1 * constants::over_c_sqr);
//    // const auto gamma2_star = particle2.gamma * com_scatter.gamma * (1.0 - com_scatter.vcm_dot_v2 * constants::over_c_sqr);
//    //
//    // auto getScatteringAngle = [&]() -> std::array<double, 2> {
//    //    const auto cos_theta = 2.0 * params.rand[3] - 1.0;
//    //    return {cos_theta, std::sqrt(1.0 - math::SQR(cos_theta))};
//    // };
//    //
//    // const auto& [p1f_star_scatter, pcoef_lab_scatter] = calculate_P1final(
//    //    com_scatter.p1_star,
//    //    com_scatter.velocity,
//    //    com_scatter.vcm_len2,
//    //    com_scatter.gamma,
//    //    params.rand[2],
//    //    getScatteringAngle
//    // );
//    //
//    // updateParticles(
//    //    particle1,
//    //    particle2,
//    //    params.mass1, params.mass2,
//    //    com_scatter.velocity,
//    //    p1f_star_scatter,
//    //    pcoef_lab_scatter,
//    //    gamma1_star, gamma2_star, com_scatter.gamma,
//    //    false, //electron_scatters,
//    //    target_ionizes
//    // );
// } // end ionizationCollision

// void coulombCollision(const bool has_coulomb, const auto& params, const auto& spec, const auto& coulomb)
// {
//    if (!has_coulomb) { return; } // no coulomb collisions
//
//    auto& particle1 = params.particle1;
//    auto& particle2 = params.particle2;
//    const auto& cspec = spec.coulomb;
//
//    const auto scatter_p1 = params.rand[0] * params.weight1 < params.weight2;
//    const auto scatter_p2 = params.rand[0] * params.weight2 < params.weight1;
//    const auto dv_length = (particle1.velocity - particle2.velocity).length();
//
//    if (dv_length == 0.0 or !(scatter_p1 or scatter_p2)) {
//       return;
//    } // parallel particles or self particles can cause divide by zero errors
//
//    const auto com = getCOMData(particle1.velocity, particle2.velocity, particle1.gamma, particle2.gamma, params.mass1, params.mass2);
//
//    const auto p1_star_len = com.p1_star.length();
//    const auto gamma1_star = particle1.gamma * com.gamma1;
//    const auto gamma2_star = particle2.gamma * com.gamma2;
//
//    const auto gamma_coef1 = com.gamma / (com.mg1 + com.mg2);
//    const auto gamma_coef2 = 1.0 + constants::c_sqr * (gamma1_star * params.mass1) * (gamma2_star * params.mass2) / math::SQR(p1_star_len);
//
//    double coulomb_log_pairwise = cspec.coulomb_log;
//    if (coulomb_log_pairwise <= 0.0) {
//       const auto l_deBroglie = 0.5 * constants::h / p1_star_len;
//       const auto b0 = gamma_coef1 * gamma_coef2 * coulomb.coef2;
//       const auto bmin2 = std::max(math::SQR(b0), math::SQR(l_deBroglie));
//       coulomb_log_pairwise = std::max(2.0, 0.5 * std::log(1.0 + coulomb.bmax2 / bmin2));
//    }
//
//    const auto s12_calc = cspec.rate_multiplier * params.max_weight * coulomb.coef1 * params.scatter_coef * gamma_coef1 * math::SQR(gamma_coef2) * coulomb_log_pairwise * p1_star_len / (particle1.gamma * particle2.gamma * constants::c_sqr);
//    const auto s12_max = coulomb.scatter_lowt * cspec.rate_multiplier * params.max_weight * dv_length;
//    const auto s12 = std::min(s12_max, s12_calc);
//
//    auto getScatteringAngles = [&]() -> std::array<double, 2> {
//       double cos_theta;
//       if (s12 >= 4.0) {
//          cos_theta = 2.0 * params.rand[1] - 1.0;
//       }
//       else {
//          const auto alpha = 0.37 * s12 - 0.005 * math::SQR(s12) - 0.0064 * math::CUBE(s12);
//          const auto sin2x2 = params.rand[1] * alpha / std::sqrt(1.0 + params.rand[1] * (math::SQR(alpha) - 1.0));
//          cos_theta = 1.0 - (2.0 * sin2x2);
//       }
//       return {cos_theta, std::sqrt(1.0 - math::SQR(cos_theta))};
//    };
//
//    const auto& [p1f_star, pcoef_lab] = calculate_P1final(
//       com.p1_star,
//       com.velocity,
//       com.vcm_len2,
//       com.gamma,
//       params.rand[2],
//       getScatteringAngles
//    );
//
//    updateParticles(
//       particle1,
//       particle2,
//       params.mass1, params.mass2,
//       com.velocity,
//       p1f_star,
//       pcoef_lab,
//       gamma1_star, gamma2_star, com.gamma,
//       scatter_p1,
//       scatter_p2
//    );
// } // end coulombCollision

} // end namespace tf::collisions

#endif //GAME_ENGINE_BINARY_CHANNELS_HPP