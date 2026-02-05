#ifndef EM_CONSTANTS_HPP
#define EM_CONSTANTS_HPP

namespace tf::constants
{
inline constexpr auto         pi{3.141592653589793};
inline constexpr auto        pi2{2.0 * pi};
inline constexpr auto          c{299792458.0};        // m/s
inline constexpr auto      c_sqr{c * c};              // (m/s)^2
inline constexpr auto over_c_sqr{1.0 / c_sqr};        // (s/m)^2

inline constexpr auto       eps0{8.854187812813e-12}; // F/m
inline constexpr auto        mu0{1.2566370621219e-6}; // H/m
inline constexpr auto       eta0{376.73031366686992}; // Ohms

inline constexpr auto        q_e{1.602176634e-19};    // electron charge, coulombs
inline constexpr auto        m_e{9.109383714e-31};    // electron mass, kg
inline constexpr auto        m_p{1.672621925e-27};    // proton mass, kg
inline constexpr auto      m_e_c{m_e * c};            // m_e * c
inline constexpr auto  m_e_c_sqr{m_e * c_sqr};        // m_e * c^2, J

inline constexpr auto         kB{1.380649e-23};       // J/K
inline constexpr auto          h{6.62607015e-34};     // Planck's Constant, J/s
} // end namespace tf::constants

#endif //EM_CONSTANTS_HPP
