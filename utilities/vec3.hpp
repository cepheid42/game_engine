#ifndef TRIFORCE_VECTOR_H
#define TRIFORCE_VECTOR_H

#include <cmath>
#include <initializer_list>
#include <sstream>

// ===== Vector Types =====
// ========================
namespace tf
{
  // ----- Detail -----
  namespace detail
  {
    template<typename T, size_t N>
    struct VectorBase {
      using value_t = T;
      // using dimension_t = tags::Dimension<N>;

      VectorBase() = default;
      VectorBase(std::initializer_list<T[N]> init_) : data(init_) {}
      
      T &operator[](int i) { return data[i]; }
      const T &operator[](int i) const { return data[i]; }
      
      auto length_squared() {
        auto sum = T(0.0);
        for (int i = 0; i < N; i++) { sum += data[i] * data[i]; }
        return sum;
      }
      T length() const { return std::sqrt(length_squared()); }

      T data[N];
      //
    };// end struct tf::detail::VectorBase
    //
  } // end namespace tf::detail

  // ===== vec3 =====
  // ================
  template<typename T>
  struct vec3 : detail::VectorBase<T, 3> {
    using value_t = typename detail::VectorBase<T, 3>::value_t;
    // using dimension_t = typename detail::VectorBase<T, 3>::dimension_t;
    
    vec3() = default;
    vec3(std::initializer_list<value_t[3]> init_) : detail::VectorBase<T, 3>(init_) {}
    vec3(value_t e0, value_t e1, value_t e2) : vec3({e0, e1, e2}) {}
    
    // Unary Negation
    vec3 operator-() const { return {-(*this)[0], -(*this)[1], -(*this)[2]}; }
    
    vec3 &operator+=(const vec3 &v) {
      (*this)[0] += v[0];
      (*this)[1] += v[1];
      (*this)[2] += v[2];
      return *this;
    }
    
    vec3 &operator-=(const vec3 &v) {
      (*this)[0] -= v[0];
      (*this)[1] -= v[1];
      (*this)[2] -= v[2];
      return *this;
    }
    
    vec3 &operator*=(const T s) {
      (*this)[0] *= s;
      (*this)[1] *= s;
      (*this)[2] *= s;
      return *this;
    }
    
    vec3 &operator/=(const T s) {
      (*this)[0] /= s;
      (*this)[1] /= s;
      (*this)[2] /= s;
      return *this;
    }
    
    bool operator==(const vec3 &v) const { return ((*this)[0] == v[0] && (*this)[1] == v[1] && (*this)[2] == v[2]); }
    bool operator!=(const vec3 &v) const { return !((*this) == v); }
  };// end struct tf::vec3
} // end namespace tf

// ===== vec3 Operators =====
// ==========================
template<typename T>
tf::vec3<T> operator*(T s, const tf::vec3<T>& u)  {
  return {s * u[0], s * u[1], s * u[2]};
}

template<typename T>
tf::vec3<T> operator*(const tf::vec3<T>& u, T s)  {
  return s * u;
}

template<typename T>
tf::vec3<T> operator/(const tf::vec3<T>& u, T s) {
  return (T(1) / s) * u;
}

template<typename T>
tf::vec3<T> operator+(const tf::vec3<T>& u, const tf::vec3<T>& v) {
  return {u[0] + v[0], u[1] + v[1], u[2] + v[2]};
}

template<typename T>
tf::vec3<T> operator-(const tf::vec3<T>& u, const tf::vec3<T>& v) {
  return {u[0] - v[0], u[1] - v[1], u[2] - v[2]};
}

template<typename T>
tf::vec3<T> operator*(const tf::vec3<T>& u, const tf::vec3<T>& v) {
  return {u[0] * v[0], u[1] * v[1], u[2] * v[2]};
}

template<typename T>
tf::vec3<T> operator/(const tf::vec3<T>& u, const tf::vec3<T>& v) {
  return {u[0] / v[0], u[1] / v[1], u[2] / v[2]};
}

template<typename T>
tf::vec3<T> unit_vector(const tf::vec3<T>& u){
  return u / u.length();
}

template<typename T>
T dot(const tf::vec3<T>& u, const tf::vec3<T>& v) {
  // Performs u @ v
  return (u[0] * v[0]) + (u[1] * v[1]) + (u[2] * v[2]);
}

template<typename T>
tf::vec3<T> cross(const tf::vec3<T>& u, const tf::vec3<T>& v) {
  // Performs u x v
  return {u[1] * v[2] - u[2] * v[1],
          u[2] * v[0] - u[0] * v[2],
          u[0] * v[1] - u[1] * v[0]};
}

///*
// *  Cartesian:   (x, y, z)
// *  Cylindrical: (r, phi, z)      phi angle around Z axis.
// *  Spherical:   (r, theta, phi)  theta angle from +Z axis.
// */
//
///* To Cartesian Coordinates */
//template<typename T>
//tf::vec3<T> cylindrical_to_cartesian(const tf::vec3<T>& v) {
//  auto x = v[0] * cos(v[1]);
//  auto y = v[0] * sin(v[1]);
//  return {x, y, v[2]};
//}
//
//template<typename T>
//tf::vec3<T> spherical_to_cartesian(const tf::vec3<T>& v) {
//  auto x = v[0] * sin(v[2]) * cos(v[1]);
//  auto y = v[0] * sin(v[2]) * sin(v[1]);
//  auto z = v[0] * cos(v[2]);
//  return {x, y, z};
//}
//
///* To Cylindrical Coordinates */
//template<typename T>
//tf::vec3<T> cartesian_to_cylindrical(const tf::vec3<T>& v) {
//  auto r = sqrt((v[0] * v[0]) + (v[1] * v[1]));
//  auto phi = atan2(v[1], v[2]);
//  return {r, phi, v[2]};
//}
//
//template<typename T>
//tf::vec3<T> spherical_to_cylindrical(const tf::vec3<T>& v) {
//  auto r = v[0] * sin(v[2]);
//  auto z = v[0] * cos(v[2]);
//  return {r, v[1], z};
//}
//
///* To Spherical Coordinates */
//template<typename T>
//tf::vec3<T> cartesian_to_spherical(const tf::vec3<T>& v) {
//  auto rho = v.length();
//  auto theta = acos(v[2] / rho);
//  auto phi = acos(v[0] / sqrt((v[0] * v[0]) + (v[1] * v[1])));
//  return {rho, theta, phi};
//}
//
//template<typename T>
//tf::vec3<T> cylindrical_to_spherical(const tf::vec3<T>& v) {
//  auto rho = sqrt((v[0] * v[0]) * (v[2] * v[2]));
//  auto theta = atan2(v[0], v[2]);
//  return {rho, theta, v[1]};
//}

template<typename T>
tf::vec3<T> rotateAround_X(const tf::vec3<T>& v, T theta) {
  // rotates tf::vec3 v around x-axis by theta radians
  auto yp = v[1] * cos(theta) - v[2] * sin(theta);
  auto zp = v[1] * sin(theta) + v[2] * cos(theta);
  return {v[0], yp, zp};
}

// rotates tf::vec3 v around y-axis by theta radians
template<typename T>
tf::vec3<T> rotateAround_Y(const tf::vec3<T>& v, T theta) {
  // rotates tf::vec3 v around y-axis by theta radians
  auto xp = v[0] * cos(theta) + v[2] * sin(theta);
  auto zp = -v[0] * sin(theta) + v[2] * cos(theta);
  return {xp, v[1], zp};
}

template<typename T>
tf::vec3<T> rotateAround_Z(const tf::vec3<T>& v, T theta) {
  // rotates tf::vec3 v around z-axis by theta radians
  auto xp = v[0] * cos(theta) - v[1] * sin(theta);
  auto yp = v[0] * sin(theta) + v[1] * cos(theta);
  return {xp, yp, v[2]};
}

template<typename T>
tf::vec3<T> rotate(const tf::vec3<T>& v, T gamma, T beta, T alpha) {
  auto v1 = rotateAround_X(v, gamma);
  auto v2 = rotateAround_Y(v1, beta);
  auto v3 = rotateAround_Z(v2, alpha);
  return v3;
}

template<typename T>
std::istringstream& operator>>(std::istringstream& in, tf::vec3<T>& v)
{
  for (std::size_t i = 0; i < 3; i++) {
    in >> v[i];
  }

  return in;
}
#endif //TRIFORCE_VECTOR_H
