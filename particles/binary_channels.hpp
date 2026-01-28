#ifndef GAME_ENGINE_BINARY_CHANNELS_HPP
#define GAME_ENGINE_BINARY_CHANNELS_HPP

#include "constants.hpp"
#include "math_utils.hpp"
#include "particles.hpp"
#include "vec3.hpp"

#include <algorithm>
#include <cmath>
#include <span>
#include <tuple>

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
   std::span<double> rand;
};
            
struct CoulombData {
   double coef1;
   double coef2;
   double bmax2;
   double scatter_lowt;
};

struct COMData {
   vec3<double> p1;
   vec3<double> vcm;
   double gamma;
   double gamma1;
   double gamma2;
   double gamma_rel;
};

auto calculate_P1f_cm(const auto& com, const auto U, const auto scatteringAngles)
-> std::tuple<vec3<double>, double>
{
   const auto& p1 = com.p1;
   const auto p1_len = p1.length();

   const auto [cos_theta, sin_theta] = scatteringAngles;
   const auto phi = 2.0 * constants::pi * U;
   const auto dv0 = sin_theta * std::cos(phi);
   const auto dv1 = sin_theta * std::sin(phi);

   // Perez (2012) eq 12
   const auto p1_perp = std::sqrt(math::SQR(p1[0]) + math::SQR(p1[1]));
   const vec3 p1f_cm {
       dv0 * (p1[0] * p1[2] / p1_perp) - dv1 * (p1[1] * p1_len / p1_perp) + cos_theta * p1[0],
       dv0 * (p1[1] * p1[2] / p1_perp) + dv1 * (p1[0] * p1_len / p1_perp) + cos_theta * p1[1],
      -dv0 * p1_perp + cos_theta * p1[2]
   };

   const auto pcoef_lab = (com.gamma - 1.0) * dot(com.vcm, p1f_cm) / com.vcm.length_squared();
   return {p1f_cm, pcoef_lab};
}

void updateParticle(auto& p, const auto m, const auto& pf, const auto pcoef, const auto& vcm, const auto gamma, const auto gamma_cm) {
   // Perez (2012) eq 13
   const auto p1f = pf + vcm * (pcoef + m * gamma * gamma_cm);
   const auto gamma1f = particles::calculateGammaP(p1f, m);
   p.velocity = p1f / (m * gamma1f);
   p.gamma = gamma1f;
}

auto getCOMData(const auto& particle1, const auto& particle2, const auto m1, const auto m2)
-> COMData
{
   const auto& v1 = particle1.velocity;
   const auto& v2 = particle2.velocity;
   const auto& gamma1 = particle1.gamma;
   const auto& gamma2 = particle2.gamma;

   const auto p1 = gamma1 * m1 * v1;
   const auto p2 = gamma2 * m2 * v2;

   const auto vcm = (p1 + p2) / (m1 * gamma1 + m2 * gamma2);
   const auto vcm2 = vcm.length_squared();
   const auto gamma_cm = 1.0 / std::sqrt(1.0 - vcm2 * constants::over_c_sqr);

   const auto u1 = p1 / m1;
   const auto u2 = p2 / m2;

   const auto vcm_dot_u1 = dot(vcm, u1);
   const auto vcm_dot_u2 = dot(vcm, u2);

   const auto gamma1_cm = (gamma1 - vcm_dot_u1 * constants::over_c_sqr) * gamma_cm;
   const auto gamma2_cm = (gamma2 - vcm_dot_u2 * constants::over_c_sqr) * gamma_cm;

   const auto u1_cm = u1 + ((gamma_cm - 1.0) / vcm2 * vcm_dot_u1 - gamma_cm * gamma1) * vcm;
   const auto u2_cm = u2 + ((gamma_cm - 1.0) / vcm2 * vcm_dot_u2 - gamma_cm * gamma2) * vcm;

   const auto dv2_cm = (u1_cm / gamma1_cm - u2_cm / gamma2_cm).length_squared();
   const auto vrel2_cm = dv2_cm / math::SQR(1.0 - dv2_cm * constants::over_c_sqr);

   const auto gamma_rel = 1.0 / std::sqrt(1.0 - vrel2_cm * constants::over_c_sqr);

   const auto p1_star = p1 + m1 * gamma1 * vcm * ((gamma_cm - 1.0) * vcm_dot_u1 / vcm2 - gamma_cm);

   return {p1_star, vcm, gamma_cm, gamma1_cm, gamma2_cm, gamma_rel};
} // end getCOMData()

// =====================================================================================================
// =============== Coulomb Collisions ==================================================================
void coulombCollision(const bool has_coulomb, const auto& params, const auto& spec, const auto& cell_data)
{
   if (!has_coulomb) { return; } // no cell_data collisions

   auto& particle1 = params.particle1;
   auto& particle2 = params.particle2;
   const auto& v1 = particle1.velocity;
   const auto& v2 = particle2.velocity;
   const auto& gamma1 = particle1.gamma;
   const auto& gamma2 = particle2.gamma;
   const auto& w1 = particle1.weight;
   const auto& w2 = particle2.weight;
   const auto& m1 = params.mass1;
   const auto& m2 = params.mass2;
   const auto& cspec = spec.coulomb;

   const auto scatter_p1 = params.rand[0] * w1 < w2;
   const auto scatter_p2 = params.rand[0] * w2 < w1;
   const auto dv_length = (v1 - v2).length();

   // parallel particles or self particles can cause divide by zero errors
   if (!(scatter_p1 or scatter_p2) or dv_length == 0.0) { return; }

   const auto com = getCOMData(particle1, particle2, m1, m2);

   const auto p1_star_len = com.p1.length();
   const auto gamma_coef1 = com.gamma / (m1 * gamma1 + m2 * gamma2);
   const auto gamma_coef2 = 1.0 + constants::c_sqr * (com.gamma1 * m1) * (com.gamma2 * m2) / math::SQR(p1_star_len);

   double coulomb_log_pairwise{cspec.coulomb_log};
   if (coulomb_log_pairwise <= 0.0) {
      const auto l_deBroglie = 0.5 * constants::h / p1_star_len;
      const auto b0 = gamma_coef1 * gamma_coef2 * cell_data.coef2;
      const auto bmin2 = std::max(math::SQR(b0), math::SQR(l_deBroglie));
      coulomb_log_pairwise = std::max(2.0, 0.5 * std::log(1.0 + cell_data.bmax2 / bmin2));
   }

   const auto s12_calc = cspec.rate_multiplier * params.max_weight * cell_data.coef1 * params.scatter_coef * gamma_coef1 * math::SQR(gamma_coef2) * coulomb_log_pairwise * p1_star_len * constants::over_c_sqr / (gamma1 * gamma2);
   const auto s12_max = cell_data.scatter_lowt * cspec.rate_multiplier * params.max_weight * dv_length;
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

   const auto& [p1f, pcoef] = calculate_P1f_cm(com, params.rand[2], getScatteringAngles());

   if (scatter_p1) {
      updateParticle(particle1, m1, p1f, pcoef, com.vcm, com.gamma1, com.gamma);
   }

   if (scatter_p2) {
      updateParticle(particle2, m2, -p1f, -pcoef, com.vcm, com.gamma2, com.gamma);
   }
} // end coulombCollision()

// =====================================================================================================
// =============== Ionization Reactions ================================================================
void ionizationCollision(
   const bool has_ionization,
   const auto& params,
   const auto& spec,
   auto& buffers,
   const auto ndups,
   const auto& cs_table
) {
   if (!has_ionization) { return; } // no ionization collisions

   auto& particle1 = params.particle1;
   auto& particle2 = params.particle2;
   const auto& v1 = particle1.velocity;
   const auto& v2 = particle2.velocity;
   const auto& gamma1 = particle1.gamma;
   const auto& gamma2 = particle2.gamma;
   const auto& w1 = particle1.weight; // Do not use weight / dup here;
   const auto& w2 = particle2.weight;
   const auto& m1 = params.mass1;
   const auto& m2 = params.mass2;
   const auto& ionization = spec.ionization;

   const auto wi_eff = w2 / ionization.rejection_multiplier;
   const auto dv_len = (v1 - v2).length();

   auto target_ionizes = w1 >= wi_eff or (w1 / wi_eff >= params.rand[2]);
   auto electron_scatters = false; //w1 < wi_eff or (wi_eff / w1 >= params.rand[2]); // todo: this is hard coded to false for ionization test

   // parallel particles can cause divide by zero errors
   if (!(target_ionizes or electron_scatters) or dv_len == 0.0) { return; }

   const auto mu = m1 * m2 / (m1 + m2); // todo: move out of here

   const auto com = getCOMData(particle1, particle2, m1, m2);
   const auto energy_com_eV = (com.gamma_rel - 1.0) * mu * constants::c_sqr / constants::q_e;

   // Skip if COM energy is below ionization energy threshold
   if (energy_com_eV <= ionization.ionization_energy) { return; }

   auto cross_section{ionization.constant_cross_section};
   if (cross_section == 0.0) {
      if (cs_table.is_outofbounds(energy_com_eV)) { return; }
      cross_section = cs_table.lerp(energy_com_eV);
   }

   const auto probability_coef = cross_section * params.max_weight * params.scatter_coef * dv_len * ionization.rate_multiplier;
   auto prod_mult = ionization.production_multiplier;
   auto scatter_probability = probability_coef * prod_mult;

   while (scatter_probability > spec.probability_search_area and prod_mult > 1.0) {
      prod_mult /= 2.0;
      scatter_probability = probability_coef * prod_mult;
   }

   // Skip if scatter probability is below random number
   if (params.rand[0] > scatter_probability) { return; }

   const auto E_scatter = params.rand[1] * (energy_com_eV - ionization.ionization_energy);
   const auto p_scatter = v1 * gamma1 * m1 * std::sqrt(E_scatter / energy_com_eV); // todo: this can be simplified a bit, dont need p for anything but getting v
   const auto gamma_scatter = particles::calculateGammaP(p_scatter, m1);
   const auto v_scatter = p_scatter / (gamma_scatter * m1);

   if (target_ionizes) {
      const auto product_weight = std::min(w2, wi_eff * ndups / ionization.production_multiplier);
      const auto E_ejected = energy_com_eV - ionization.ionization_energy - E_scatter;
      const auto p_ejected = v1 * gamma1 * m1 * std::sqrt(E_ejected / energy_com_eV); // todo: this can be simplified a bit, dont need p for anything but getting v
      const auto gamma_ejected = particles::calculateGammaP(p_ejected, m1);
      const auto v_ejected = p_ejected / (gamma_ejected * m1);

      // todo: is this enough to keep OpenMP threads from stepping on each other?
      #pragma omp critical
      {
         buffers.g1_products.emplace_back(
            v_ejected,
            gamma_ejected,
            particle2.location,
            particle2.location,
            product_weight
         );

         buffers.g2_products.emplace_back(
            particle2.velocity,
            particle2.gamma,
            particle2.location,
            particle2.location,
            product_weight
         );
      }

      particle2.weight -= product_weight;
   } // end target_ionizes

   particles::Particle e_scattered{
      v_scatter,
      gamma_scatter,
      {}, {}, -1.0
   };

   auto getScatteringAngle = [&]() -> std::array<double, 2> {
      const auto cos_theta = 2.0 * params.rand[4] - 1.0;
      return {cos_theta, std::sqrt(1.0 - math::SQR(cos_theta))};
   };

   const auto com_scatter = getCOMData(e_scattered, particle2, m1, m2);
   const auto& [p1f, pcoef] = calculate_P1f_cm(com_scatter, params.rand[3], getScatteringAngle());

   if (electron_scatters) {
      updateParticle(particle1, m1, p1f, pcoef, com_scatter.vcm, com_scatter.gamma1, com_scatter.gamma);
   }

   if (target_ionizes) {
      updateParticle(particle2, m2, -p1f, -pcoef, com_scatter.vcm, com_scatter.gamma2, com_scatter.gamma);
   }
} // end ionizationCollision()

// =====================================================================================================
// =============== Fusion Reactions ====================================================================
void fusionCollision(
   const bool has_fusion,
   const auto& params,
   const auto& spec,
   auto& buffers,
   const auto& cs_table,
   const auto prod_m1, const auto prod_m2
) {
   if (!has_fusion) { return; } // no fusion collisions

   auto& particle1 = params.particle1;
   auto& particle2 = params.particle2;
   const auto& v1 = particle1.velocity;
   const auto& v2 = particle2.velocity;
   const auto& w1 = particle1.weight;
   const auto& w2 = particle2.weight;
   // const auto& w1 = params.weight1;
   // const auto& w2 = params.weight2; // todo: does this affect DT fusion? It doesn't affect DD.
   const auto& m1 = params.mass1;
   const auto& m2 = params.mass2;
   const auto& gamma1 = particle1.gamma;
   const auto& gamma2 = particle2.gamma;
   const auto& fusion = spec.fusion;

   const auto dv = v1 - v2;
   const auto dv_len = dv.length();

   const auto mu = m1 * m2 / (m1 + m2); // todo: move out of here

   const auto p1 = gamma1 * m1 * v1;
   const auto p2 = gamma2 * m2 * v2;

   const auto vcm = (p1 + p2) / (m1 * gamma1 + m2 * gamma2);
   const auto vcm2 = vcm.length_squared();
   const auto gamma_cm = 1.0 / std::sqrt(1.0 - vcm2 * constants::over_c_sqr);

   const auto u1 = p1 / m1;
   const auto u2 = p2 / m2;

   const auto vcm_dot_u1 = dot(vcm, u1);
   const auto vcm_dot_u2 = dot(vcm, u2);

   const auto gamma1_cm = (gamma1 - vcm_dot_u1 * constants::over_c_sqr) * gamma_cm;
   const auto gamma2_cm = (gamma2 - vcm_dot_u2 * constants::over_c_sqr) * gamma_cm;

   const auto u1_cm = u1 + vcm * ((gamma_cm - 1.0) / vcm2 * vcm_dot_u1 - gamma_cm * gamma1);
   const auto u2_cm = u2 + vcm * ((gamma_cm - 1.0) / vcm2 * vcm_dot_u2 - gamma_cm * gamma2);

   const auto dv2_cm = (u1_cm / gamma1_cm - u2_cm / gamma2_cm).length_squared();
   const auto vrel2_cm = dv2_cm / math::SQR(1.0 - dv2_cm * constants::over_c_sqr);

   const auto gamma_rel = 1.0 / std::sqrt(1.0 - vrel2_cm * constants::over_c_sqr);
   const auto energy_com_eV = (gamma_rel - 1.0) * mu * constants::c_sqr / constants::q_e;

   if (cs_table.is_outofbounds(energy_com_eV) or dv_len == 0.0) { return; } // Skip if energy is below minimum of table

   auto cross_section{fusion.constant_cross_section};
   if (cross_section == 0.0) {
      cross_section = cs_table.lerp(energy_com_eV);
   }

   const auto probability_coef = cross_section * params.max_weight * params.scatter_coef * dv_len * fusion.rate_multiplier;
   auto prod_mult = fusion.production_multiplier;
   auto scatter_probability = probability_coef * prod_mult;

   while (scatter_probability > spec.probability_search_area and prod_mult > 1.0) {
      prod_mult /= 10.0;
      scatter_probability = probability_coef * prod_mult;
   }

   if (params.rand[0] > scatter_probability) { return; }

   const auto dv_perp = std::sqrt(math::SQR(dv[0]) + math::SQR(dv[1]));
   const auto cos_theta = dv[2] / dv_len;
   const auto sin_theta = dv_perp / dv_len;

   const auto cos_phi = dv[0] / dv_perp;
   const auto sin_phi = dv[1] / dv_perp;

   const vec3 vp{
      cos_phi * cos_theta * dv[0] + sin_phi * cos_theta * dv[1] - sin_theta * dv[2],
      -sin_phi * dv[0] + cos_phi * dv[1],
      cos_phi * sin_theta * dv[0] + sin_phi * sin_theta * dv[1] + cos_theta * dv[2]
   };

   const auto v_cm = vp * mu / m2;

   const auto phi_scatter = 2.0 * constants::pi * params.rand[1];
   const auto sin_phi_scatter = std::sin(phi_scatter);
   const auto cos_phi_scatter = std::cos(phi_scatter);

   auto getScatteringAngle = [&]() -> std::array<double, 2> {
      const auto ctheta = 2.0 * params.rand[4] - 1.0;
      return {ctheta, std::sqrt(1.0 - math::SQR(ctheta))};
   };

   const auto& [cos_theta_scatter, sin_theta_scatter] = getScatteringAngle();
   const vec3 rotation_vec {
      sin_theta_scatter * cos_phi_scatter,
      sin_theta_scatter * sin_phi_scatter,
      cos_theta_scatter
   };

   const auto mass_ratio = prod_m1 / prod_m2;
   const auto v3_rot_com_len = std::sqrt(2.0 / (mass_ratio * (prod_m1 + prod_m2)) * constants::q_e * (energy_com_eV + fusion.energy_gain));

   const auto v3_rot_cm = v3_rot_com_len * rotation_vec;
   const auto v4_rot_cm = -1.0 * mass_ratio * v3_rot_cm; // Mike says the sqrt in the paper is incorrect

   const auto v3_rot = v3_rot_cm + v_cm;
   const auto v4_rot = v4_rot_cm + v_cm;

   const vec3 u3{
      cos_phi * cos_theta * v3_rot[0] - sin_phi * v3_rot[1] + cos_phi * sin_theta * v3_rot[2],
      sin_phi * cos_theta * v3_rot[0] + cos_phi * v3_rot[1] + sin_phi * sin_theta * v3_rot[2],
      -sin_theta * v3_rot[0] + cos_theta * v3_rot[2]
   };

   const vec3 u4{
      cos_phi * cos_theta * v4_rot[0] - sin_phi * v4_rot[1] + cos_phi * sin_theta * v4_rot[2],
      sin_phi * cos_theta * v4_rot[0] + cos_phi * v4_rot[1] + sin_phi * sin_theta * v4_rot[2],
      -sin_theta * v4_rot[0] + cos_theta * v4_rot[2]
   };

   const auto v_prod1 = u3 + v2;
   const auto v_prod2 = u4 + v2;

   const auto gamma_prod1 = particles::calculateGammaV(v_prod1);
   const auto gamma_prod2 = particles::calculateGammaV(v_prod2);

   const auto product_weight = std::min(w1, w2) / prod_mult;

   #pragma omp critical
   {
      // Half of each product is placed at each input particle position, for 4 total particles
      buffers.g1_products.emplace_back(
         v_prod1,
         gamma_prod1,
         particle1.location,
         particle1.location,
         product_weight / 2.0
      );

      buffers.g1_products.emplace_back(
         v_prod1,
         gamma_prod1,
         particle2.location,
         particle2.location,
         product_weight / 2.0
      );

      buffers.g2_products.emplace_back(
         v_prod2,
         gamma_prod2,
         particle1.location,
         particle1.location,
         product_weight / 2.0
      );

      buffers.g2_products.emplace_back(
         v_prod2,
         gamma_prod2,
         particle2.location,
         particle2.location,
         product_weight / 2.0
      );
   }

   particle1.weight -= product_weight;
   particle2.weight -= product_weight;
} // end fusionCollision()


// void bremsstrahlungCollision(
//    bool has_brem,
//    const auto& params,
//    const auto& spec,
//    auto& buffers,
//    const auto& cs_table,
//    const auto prod_m1, const auto prod_m2
// )
// {
//    if (!has_brem) { return; }
//
//    auto& particle1 = params.particle1;
//    auto& particle2 = params.particle2;
//    const auto& v1 = particle1.velocity;
//    const auto& v2 = particle2.velocity;
//    // const auto& w1 = particle1.weight;
//    // const auto& w2 = particle2.weight;
//    const auto& w1 = params.weight1;
//    const auto& w2 = params.weight2; // todo: does this affect DT fusion? It doesn't affect DD.
//    const auto& m1 = params.mass1;
//    const auto& m2 = params.mass2;
//    const auto& gamma1 = particle1.gamma;
//    const auto& gamma2 = particle2.gamma;
//    const auto& brem = spec.brem;
//
//    const auto beta_ve = v1 / constants::c;
//    const auto beta_vi = v2 / constants::c;
//    const auto beta_i2 = beta_vi.length_squared();
//
//    const auto beta_ei = dot(beta_ve, beta_vi);
//    const auto gamma_ei = gamma1 * gamma2;
//
//    // Calculation is done in ion frame, so compute lorentz transform of electron velocity
//    const auto coef = (gamma2 - 1.0) * beta_ei / beta_i2 - gamma2;
//    const auto p_ep = gamma1 * m1 * v1 + coef * constants::m_e * constants::c_sqr * gamma1 * beta_ve;
//    const auto gamma_ep = gamma_ei * (1.0 - beta_ei);
//    const auto v_ep = p_ep.length() / (constants::m_e * gamma_ep);
//
//    // Electron energy in ion frame (eV)
//    const auto electron_energy_eV = (gamma_ep - 1.0) * constants::m_e * constants::c_sqr / constants::q_e;
//
//    // todo: differential cross section
//    // Total cross-section and photon energy from differential cross-section
//    const auto cross_section_m2 = cs_table.lerp(electron_energy_eV);
//
//    const auto probability_coef = cross_section_m2 * params.max_weight * params.scatter_coef * brem.rate_multiplier * v_ep * gamma_ep / gamma_ei;
//
//    auto prod_mult = brem.production_multiplier;
//    auto scatter_probability = probability_coef * prod_mult;
//
//    // Method requires P<1 so reduced the multiplier until that is satisfied
//    while (scatter_probability > 1.0 and prod_mult > 1.0) {
//       prod_mult /= 10.0;
//       scatter_probability = probability_coef * prod_mult;
//    }
//
//    // Production multiplier only applies to likelihood of photon production
//    // and reduces weight accordingly. Use the probability computed without
//    // the multiplier to determine if the electron loses energy
//    auto make_photon = params.rand[0] < scatter_probability;
//    auto remove_electron_energy = params.rand[0] < probability_coef and params.rand[0] <= w2 / params.max_weight and reduce_electron_energy;
//
//    if (!(make_photon or remove_electron_energy)) { return; }
//
//    // Sample photon energy from cumulative differential cross-section
//    // start by finding the closest electron energy in table.
//
//    // todo: here be dragons
//
// } // end bremsstrahlungCollision()

} // end namespace tf::collisions

#endif //GAME_ENGINE_BINARY_CHANNELS_HPP