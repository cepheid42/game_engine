#ifndef EM_CONSTANTS_HPP
#define EM_CONSTANTS_HPP

namespace tf::constants
{
  // Pi values
  template<typename T> constexpr T pi = static_cast<T>(3.141592653589793);
  template<typename T> constexpr T pi2 = static_cast<T>(2.0 * pi<T>);

  // Useful
  template<typename T> constexpr T kB = static_cast<T>(1.380649e-23); // J/K
  // template<typename T> constexpr T root2 = 1.414213562373095);
  // template<typename T> constexpr T root3 = 1.732050807568877);

  // Physical
  template<typename T> constexpr T c = static_cast<T>(299792458.0);  // m/s
  template<typename T> constexpr T c_sqr = static_cast<T>(c<T> * c<T>);    // (m/s)^2
  template<typename T> constexpr T over_c_sqr = static_cast<T>(1.0 / c_sqr<T>);    // (s/m)^2

  template<typename T> constexpr T eps0 = static_cast<T>(8.854187812813e-12);   // F/m
  template<typename T> constexpr T mu0 = static_cast<T>(1.2566370621219e-6);    // H/m
  template<typename T> constexpr T eta0 = static_cast<T>(376.73031366686992);   // Ohms

  template<typename T> constexpr T q_e = static_cast<T>(1.6021766e-19);        // electron charge, coulombs
  template<typename T> constexpr T m_e = static_cast<T>(9.10938370e-31);       // electron mass, kg
  template<typename T> constexpr T m_p = static_cast<T>(1.67262193e-27);       // proton mass, kg
}

#endif //EM_CONSTANTS_HPP
