#ifndef TRIFORCE_VECTOR_H
#define TRIFORCE_VECTOR_H

#include <cmath>
#include <functional>
#include <sstream>
#include <string>
#include <format>
#include <print>

// ===== Vector Types =====
// ========================
namespace tf
{
  // ===== vec3 =====
  // ================
  template<typename T>
  struct vec3 {    
    constexpr vec3() = default;
    constexpr vec3(T e0, T e1, T e2) : data{e0, e1, e2} {}
    
    T &operator[](std::size_t i) { return data[i]; }
    const T &operator[](std::size_t i) const { return data[i]; }
      
    [[nodiscard]] auto length_squared() const {
      auto sum = T(0.0);
      for (std::size_t i = 0; i < 3; i++) { sum += data[i] * data[i]; }
      return sum;
    }
    
    [[nodiscard]] T length() const { return std::sqrt(length_squared()); }

    // Unary Negation
    vec3 operator-() const { return {-data[0], -data[1], -data[2]}; }
    
    vec3& operator+=(const vec3 &v) {
      data[0] += v[0];
      data[1] += v[1];
      data[2] += v[2];
      return *this;
    }
    
    vec3& operator-=(const vec3 &v) {
      data[0] -= v[0];
      data[1] -= v[1];
      data[2] -= v[2];
      return *this;
    }
    
    vec3& operator*=(const T s) {
      data[0] *= s;
      data[1] *= s;
      data[2] *= s;
      return *this;
    }
    
    vec3& operator/=(const T s) {
      data[0] /= s;
      data[1] /= s;
      data[2] /= s;
      return *this;
    }

    vec3& operator/=(const vec3& v) {
      data[0] /= v[0];
      data[1] /= v[1];
      data[2] /= v[2];
      return *this;
    }

    template<typename U>
    vec3<U> as_type() const {
      return {static_cast<U>(data[0]),
              static_cast<U>(data[1]),
              static_cast<U>(data[2])};
    }

    bool operator==(const vec3 &v) const { return (data[0] == v[0] && data[1] == v[1] && data[2] == v[2]); }
    bool operator!=(const vec3 &v) const { return !(data == v); }

    T data[3];
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

template <typename T>
struct std::formatter<tf::vec3<T>> : std::formatter<std::string> {
  auto format(const tf::vec3<T>& p, std::format_context& ctx) const {
    return std::format_to(ctx.out(), "[{}, {}, {}]", p[0], p[1], p[2]);
  }
};

#endif //TRIFORCE_VECTOR_H
