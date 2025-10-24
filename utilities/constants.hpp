#ifndef EM_CONSTANTS_HPP
#define EM_CONSTANTS_HPP

namespace tf::constants
{
template <typename T> inline constexpr auto         pi = 3.141592653589793;
template <typename T> inline constexpr auto        pi2 = 2.0 * pi<T>;
template <typename T> inline constexpr auto          c = 299792458.0;        // m/s
template <typename T> inline constexpr auto      c_sqr = c<T> * c<T>;        // (m/s)^2
template <typename T> inline constexpr auto over_c_sqr = 1.0 / c_sqr<T>;     // (s/m)^2
template <typename T> inline constexpr auto       eps0 = 8.854187812813e-12; // F/m
template <typename T> inline constexpr auto        mu0 = 1.2566370621219e-6; // H/m
template <typename T> inline constexpr auto       eta0 = 376.73031366686992; // Ohms
template <typename T> inline constexpr auto        q_e = 1.602176634e-19;    // electron charge, coulombs
template <typename T> inline constexpr auto        m_e = 9.109383714e-31;    // electron mass, kg
template <typename T> inline constexpr auto        m_p = 1.672621925e-27;    // proton mass, kg
template <typename T> inline constexpr auto         kB = 1.380649e-23;       // J/K
template <typename T> inline constexpr auto          h = 6.62607015e-34;     // Planck's Constant, J/s
} // end namespace tf::constants
// template<typename T> inline constexpr auto root2 = 1.414213562373095;
// template<typename T> inline constexpr auto root3 = 1.732050807568877;
#endif //EM_CONSTANTS_HPP
