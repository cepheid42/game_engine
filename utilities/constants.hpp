#ifndef EM_CONSTANTS_HPP
#define EM_CONSTANTS_HPP

namespace tf::constants
{
  // Pi values
  inline constexpr double pi { 3.141592653589793 };
  inline constexpr double pi2 { 2.0 * pi };
  inline constexpr double pi4 { 4.0 * pi };
  inline constexpr double pi6 { 6.0 * pi };
  inline constexpr double pi_sqr { pi * pi };
  inline constexpr double pi_half { pi / 2.0 };
  inline constexpr double over_pi { 1.0 / pi };

  // Useful
  inline constexpr double root2 { 1.414213562373095 };
  inline constexpr double root3 { 1.732050807568877 };
  inline constexpr double one_half { 1.0 / 2.0 };
  inline constexpr double one_third { 1.0 / 3.0 };
  inline constexpr double one_quarter { 1.0 / 4.0 };
  inline constexpr double one_fifth { 1.0 / 5.0 };
  inline constexpr double one_sixth { 1.0 / 6.0 };

  // Physical
  inline constexpr double c { 299792458.0 };  // m/s
  inline constexpr double c_sqr { c * c };    // (m/s)^2
  inline constexpr double over_c_sqr { 1 / (c * c) };    // (s/m)^2

  inline constexpr double eps0 { 8.854187812813e-12 };   // F/m
  inline constexpr double mu0 { 1.2566370621219e-6 };        // H/m
  inline constexpr double eta0 { 376.73031366686992 };       // Ohms
  inline constexpr double sqrt_mueps { 3.3356409519864e-9 }; // sqrt(mu_0 * eps_0)
  inline constexpr double four_pi_eps0 { 1.11265005544787e-10 };

}

#endif //EM_CONSTANTS_HPP
