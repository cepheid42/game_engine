#ifndef EM_CONSTANTS_HPP
#define EM_CONSTANTS_HPP

namespace tf::constants
{
  // Pi values
  inline constexpr auto pi { 3.141592653589793 };
  inline constexpr auto pi2 { 2.0 * pi };
  // inline constexpr auto pi4 { 4.0 * pi };
  // inline constexpr auto pi6 { 6.0 * pi };
  // inline constexpr auto pi_sqr { pi * pi };
  // inline constexpr auto pi_half { pi / 2.0 };
  // inline constexpr auto over_pi { 1.0 / pi };

  // Useful
  inline constexpr auto root2 { 1.414213562373095 };
  inline constexpr auto root3 { 1.732050807568877 };

  // Physical
  inline constexpr auto c { 299792458.0 };  // m/s
  inline constexpr auto c_sqr { c * c };    // (m/s)^2
  inline constexpr auto over_c_sqr { 1.0 / (c * c) };    // (s/m)^2

  inline constexpr auto eps0 { 8.854187812813e-12 };   // F/m
  inline constexpr auto mu0 { 1.2566370621219e-6 };    // H/m
  inline constexpr auto eta0 { 376.73031366686992 };   // Ohms
  // inline constexpr auto sqrt_mueps { 3.3356409519864e-9 }; // sqrt(mu_0 * eps_0)
  // inline constexpr auto four_pi_eps0 { 1.11265005544787e-10 };

  inline constexpr auto q_e { 1.6021766e-19 };        // electron charge, coulombs
  inline constexpr auto m_e { 9.10938370e-31 };       // electron mass, kg

}

#endif //EM_CONSTANTS_HPP
