#ifndef EM_CONSTANTS_HPP
#define EM_CONSTANTS_HPP

namespace tf::constants
{
  // Pi values
  inline constexpr auto pi { 3.141592653589793f };
  inline constexpr auto pi2 { 2.0f * pi };
  inline constexpr auto pi4 { 4.0f * pi };
  inline constexpr auto pi6 { 6.0f * pi };
  inline constexpr auto pi_sqr { pi * pi };
  inline constexpr auto pi_half { pi / 2.0f };
  inline constexpr auto over_pi { 1.0f / pi };

  // Useful
  inline constexpr auto root2 { 1.414213562373095f };
  inline constexpr auto root3 { 1.732050807568877f };
  inline constexpr auto one_half { 1.0f / 2.0f };
  inline constexpr auto one_third { 1.0f / 3.0f };
  inline constexpr auto one_quarter { 1.0f / 4.0f };
  inline constexpr auto one_fifth { 1.0f / 5.0f };
  inline constexpr auto one_sixth { 1.0f / 6.0f };

  // Physical
  inline constexpr auto c { 299792458.0f };  // m/s
  inline constexpr auto c_sqr { c * c };    // (m/s)^2
  inline constexpr auto over_c_sqr { 1.0f / (c * c) };    // (s/m)^2

  inline constexpr auto eps0 { 8.854187812813e-12f };   // F/m
  inline constexpr auto mu0 { 1.2566370621219e-6f };        // H/m
  inline constexpr auto eta0 { 376.73031366686992f };       // Ohms
  inline constexpr auto sqrt_mueps { 3.3356409519864e-9f }; // sqrt(mu_0 * eps_0)
  inline constexpr auto four_pi_eps0 { 1.11265005544787e-10f };

}

#endif //EM_CONSTANTS_HPP
