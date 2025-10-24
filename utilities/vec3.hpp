#ifndef TRIFORCE_VECTOR_H
#define TRIFORCE_VECTOR_H

#include <cmath>
#include <utility>
#include <string>
#include <iostream>

// ===== Vector Types =====
// ========================
namespace tf {
// ===== vec3 =====
// ================
template <typename T>
struct vec3 {
   using type = T;

   T x, y, z;

   constexpr auto operator[](std::size_t i)       ->       T& { return *(reinterpret_cast<T*>(this) + i); }
   constexpr auto operator[](std::size_t i) const -> const T& { return *(reinterpret_cast<const T*>(this) + i); }

   [[nodiscard]] constexpr auto length_squared() const -> T { return x * x + y * y + z * z; }
   [[nodiscard]] constexpr auto length() const -> T { return std::hypot(x, y, z); }

   // Unary Negation
   constexpr auto operator-() const -> vec3 { return {-x, -y, -z}; }

   constexpr auto operator+=(const vec3& v) -> vec3& {
      x += v.x;
      y += v.y;
      z += v.z;
      return *this;
   }

   constexpr auto operator-=(const vec3& v) -> vec3& {
      x -= v.x;
      y -= v.y;
      z -= v.z;
      return *this;
   }

   constexpr auto operator*=(const T s) -> vec3& {
      x *= s;
      y *= s;
      z *= s;
      return *this;
   }

   constexpr auto operator/=(const T s) -> vec3& {
      x /= s;
      y /= s;
      z /= s;
      return *this;
   }

   // vec3& operator/=(const vec3& v) {
   //    this->x /= v.x;
   //    this->y /= v.y;
   //    this->z /= v.z;
   //    return *this;
   // }

   template <typename U>
   constexpr auto as_type() const -> vec3<U> {
      return vec3<U>{static_cast<U>(x),
                     static_cast<U>(y),
                     static_cast<U>(z)};
   }

   friend constexpr bool operator==(const vec3& u, const vec3& v) { return u.x == v.x and u.y == v.y and u.z == v.z; }
   friend constexpr bool operator!=(const vec3& u, const vec3& v) { return !(u == v); }

   template <std::size_t I> constexpr auto& get() & { return this->operator[](I); }
   template <std::size_t I> constexpr const auto& get() const & { return this->operator[](I); }
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
constexpr auto operator*(const T s, const tf::vec3<T>& u) -> tf::vec3<T> {
   return {s * u.x, s * u.y, s * u.z};
}

template <typename T>
constexpr auto operator*(const tf::vec3<T>& u, const T s) -> tf::vec3<T> {
   return s * u;
}

template <typename T>
constexpr auto operator/(const tf::vec3<T>& u, const T s) -> tf::vec3<T> {
   return (T(1) / s) * u;
}

template <typename T>
constexpr auto operator+(const tf::vec3<T>& u, const T s) -> tf::vec3<T> {
   return {u.x + s, u.y + s, u.z + s};
}

template <typename T>
constexpr auto operator+(const T s, const tf::vec3<T>& u) -> tf::vec3<T> {
   return u + s;
}

template <typename T>
constexpr auto operator-(const tf::vec3<T>& u, const T s) -> tf::vec3<T> {
   return {u.x - s, u.y - s, u.z - s};
}

template <typename T>
constexpr auto operator-(const T s, const tf::vec3<T>& u) -> tf::vec3<T> {
   return {s - u.x, s - u.y, s - u.z};
}

// ===== vec3-vec3 Operators =====
// =================================
template <typename T>
constexpr auto operator+(const tf::vec3<T>& u, const tf::vec3<T>& v) -> tf::vec3<T> {
   return {u.x + v.x, u.y + v.y, u.z + v.z};
}

template <typename T>
constexpr auto operator-(const tf::vec3<T>& u, const tf::vec3<T>& v) -> tf::vec3<T> {
   return {u.x - v.x, u.y - v.y, u.z - v.z};
}

template <typename T>
constexpr auto operator*(const tf::vec3<T>& u, const tf::vec3<T>& v) -> tf::vec3<T> {
   return {u.x * v.x, u.y * v.y, u.z * v.z};
}

template <typename T>
constexpr auto operator/(const tf::vec3<T>& u, const tf::vec3<T>& v) -> tf::vec3<T> {
   return {u.x / v.x, u.y / v.y, u.z / v.z};
}

template <typename T>
constexpr auto unit_vector(const tf::vec3<T>& u) -> tf::vec3<T> {
   return u / u.length();
}

template <typename T>
constexpr auto dot(const tf::vec3<T>& u, const tf::vec3<T>& v) -> T {
   // Performs u @ v
   return (u.x * v.x) + (u.y * v.y) + (u.z * v.z);
}

template <typename T>
constexpr auto cross(const tf::vec3<T>& u, const tf::vec3<T>& v) -> tf::vec3<T> {
   // Performs u x v
   return {u.y * v.z - u.z * v.y,
           u.z * v.x - u.x * v.z,
           u.x * v.y - u.y * v.x};
}

template<typename T>
constexpr auto is_equal(const tf::vec3<T>& u, const tf::vec3<T>& v) -> tf::vec3<bool> {
   return {u.x == v.x, u.y == v.y, u.z == v.z};
}


template <typename T>
constexpr auto operator>>(std::istringstream& in, tf::vec3<T>& v) -> std::istringstream& {
   for (std::size_t i = 0; i < 3; i++) {
      in >> v[i];
   }
   return in;
}

template <typename T>
constexpr auto operator<<(std::ostream& os, const tf::vec3<T>& v) -> std::ostream& {
   os << "[" << v.x << ", " << v.y << ", " << v.z << "]";
   return os;
}

template <typename T>
struct std::formatter<tf::vec3<T>> : std::formatter<std::string> {
   constexpr auto format(const tf::vec3<T>& p, std::format_context& ctx) const {
      return std::format_to(ctx.out(), "({}, {}, {})", p[0], p[1], p[2]);
   }
};

#endif //TRIFORCE_VECTOR_H
