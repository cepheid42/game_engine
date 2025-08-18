#ifndef TRIFORCE_VECTOR_H
#define TRIFORCE_VECTOR_H

#include <cmath>
#include <utility>
#include <string>
#include <print>
#include <iostream>

// ===== Vector Types =====
// ========================
namespace tf {
// ===== vec3 =====
// ================
template <typename T>
struct vec3 {
   using type = T;

   constexpr vec3() = default;

   constexpr vec3(T e0, T e1, T e2) : data{e0, e1, e2} {}

   constexpr T& operator[](std::size_t i) { return data[i]; }
   constexpr const T& operator[](std::size_t i) const { return data[i]; }

   [[nodiscard]] constexpr auto length_squared() const {
      return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
   }

   [[nodiscard]] constexpr T length() const { return std::hypot(data[0], data[1], data[2]); }

   // Unary Negation
   constexpr vec3 operator-() const { return {-data[0], -data[1], -data[2]}; }

   constexpr vec3& operator+=(const vec3& v) {
      data[0] += v[0];
      data[1] += v[1];
      data[2] += v[2];
      return *this;
   }

   constexpr vec3& operator-=(const vec3& v) {
      data[0] -= v[0];
      data[1] -= v[1];
      data[2] -= v[2];
      return *this;
   }

   constexpr vec3& operator*=(const T s) {
      data[0] *= s;
      data[1] *= s;
      data[2] *= s;
      return *this;
   }

   constexpr vec3& operator/=(const T s) {
      data[0] /= s;
      data[1] /= s;
      data[2] /= s;
      return *this;
   }

   // vec3& operator/=(const vec3& v) {
   //    data[0] /= v[0];
   //    data[1] /= v[1];
   //    data[2] /= v[2];
   //    return *this;
   // }

   template <typename U>
   constexpr vec3<U> as_type() const {
      return {
         static_cast<U>(data[0]),
         static_cast<U>(data[1]),
         static_cast<U>(data[2])
      };
   }

   friend constexpr bool operator==(const vec3& u, const vec3& v){ return (u[0] == v[0] && u[1] == v[1] && u[2] == v[2]); }
   friend constexpr bool operator!=(const vec3& u, const vec3& v) { return !(u == v); }

   template <std::size_t I>
   constexpr auto& get() & { return data[I]; }

   template <std::size_t I>
   constexpr const auto& get() const & { return data[I]; }

   T data[3];
}; // end struct tf::vec3
} // end namespace tf

namespace std {
template <typename T>
struct tuple_size<tf::vec3<T>> : std::integral_constant<std::size_t, 3> {};

template <typename T>
struct tuple_element<0, tf::vec3<T>> {
   using type = T;
};

template <typename T>
struct tuple_element<1, tf::vec3<T>> {
   using type = T;
};

template <typename T>
struct tuple_element<2, tf::vec3<T>> {
   using type = T;
};
} // end namespace std

// ===== vec3-scalar Operators =====
// =================================
template <typename T>
constexpr tf::vec3<T> operator*(const T s, const tf::vec3<T>& u) {
   return {s * u[0], s * u[1], s * u[2]};
}

template <typename T>
constexpr tf::vec3<T> operator*(const tf::vec3<T>& u, const T s) {
   return s * u;
}

template <typename T>
constexpr tf::vec3<T> operator/(const tf::vec3<T>& u, const T s) {
   return (T(1) / s) * u;
}

template <typename T>
constexpr tf::vec3<T> operator+(const tf::vec3<T>& u, const T s) {
   return {u[0] + s, u[1] + s, u[2] + s};
}

template <typename T>
constexpr tf::vec3<T> operator+(const T s, const tf::vec3<T>& u) {
   return u + s;
}

template <typename T>
constexpr tf::vec3<T> operator-(const tf::vec3<T>& u, const T s) {
   return {u[0] - s, u[1] - s, u[2] - s};
}

template <typename T>
constexpr tf::vec3<T> operator-(const T s, const tf::vec3<T>& u) {
   return {s - u[0], s - u[1], s - u[2]};
}

// ===== vec3-vec3 Operators =====
// =================================
template <typename T>
constexpr tf::vec3<T> operator+(const tf::vec3<T>& u, const tf::vec3<T>& v) {
   return {u[0] + v[0], u[1] + v[1], u[2] + v[2]};
}

template <typename T>
constexpr tf::vec3<T> operator-(const tf::vec3<T>& u, const tf::vec3<T>& v) {
   return {u[0] - v[0], u[1] - v[1], u[2] - v[2]};
}

template <typename T>
constexpr tf::vec3<T> operator*(const tf::vec3<T>& u, const tf::vec3<T>& v) {
   return {u[0] * v[0], u[1] * v[1], u[2] * v[2]};
}

template <typename T>
constexpr tf::vec3<T> operator/(const tf::vec3<T>& u, const tf::vec3<T>& v) {
   return {u[0] / v[0], u[1] / v[1], u[2] / v[2]};
}

template <typename T>
constexpr tf::vec3<T> unit_vector(const tf::vec3<T>& u) {
   return u / u.length();
}

template <typename T>
constexpr T dot(const tf::vec3<T>& u, const tf::vec3<T>& v) {
   // Performs u @ v
   return (u[0] * v[0]) + (u[1] * v[1]) + (u[2] * v[2]);
}

template <typename T>
constexpr tf::vec3<T> cross(const tf::vec3<T>& u, const tf::vec3<T>& v) {
   // Performs u x v
   return {
      u[1] * v[2] - u[2] * v[1],
      u[2] * v[0] - u[0] * v[2],
      u[0] * v[1] - u[1] * v[0]
   };
}

template <typename T>
constexpr std::istringstream& operator>>(std::istringstream& in, tf::vec3<T>& v) {
   for (std::size_t i = 0; i < 3; i++) {
      in >> v[i];
   }
   return in;
}

template <typename T>
constexpr std::ostream& operator<<(std::ostream& os, const tf::vec3<T>& v) {
   os << "[" << v[0] << ", " << v[1] << ", " << v[2] << "]";
   return os;
}

template <typename T>
struct std::formatter<tf::vec3<T>> : std::formatter<std::string> {
   constexpr auto format(const tf::vec3<T>& p, std::format_context& ctx) const {
      return std::format_to(ctx.out(), "({}, {}, {})", p[0], p[1], p[2]);
   }
};

#endif //TRIFORCE_VECTOR_H
