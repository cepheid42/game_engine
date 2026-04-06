#ifndef TRIFORCE_VECTOR_H
#define TRIFORCE_VECTOR_H

#include <cmath>
#include <format>
#include <iostream>
#include <string>
#include <type_traits>
#include <x86intrin.h>


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

   [[nodiscard]] constexpr auto to_uint() const -> vec3<std::size_t> {
      return vec3<std::size_t>{static_cast<std::size_t>(x),
                               static_cast<std::size_t>(y),
                               static_cast<std::size_t>(z)};
   }

   [[nodiscard]] constexpr auto to_float() const -> vec3<float> {
      return vec3<float>{static_cast<float>(x),
                         static_cast<float>(y),
                         static_cast<float>(z)};
   }

   [[nodiscard]] constexpr auto to_double() const -> vec3<double> {
      return vec3<double>{static_cast<double>(x),
                          static_cast<double>(y),
                          static_cast<double>(z)};
   }

   // friend constexpr bool operator==(const vec3& u, const vec3& v) { return u.x == v.x and u.y == v.y and u.z == v.z; }
   // friend constexpr bool operator!=(const vec3& u, const vec3& v) { return !(u == v); }

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
// ===============================
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
constexpr auto dot_product(const tf::vec3<T>& u, const tf::vec3<T>& v) -> T {
   // Performs u @ v
   return (u.x * v.x) + (u.y * v.y) + (u.z * v.z);
}

template <typename T>
constexpr auto cross_product(const tf::vec3<T>& u, const tf::vec3<T>& v) -> tf::vec3<T> {
   // Performs u x v
   return {u.y * v.z - u.z * v.y,
           u.z * v.x - u.x * v.z,
           u.x * v.y - u.y * v.x};
}

// template<typename T>
// requires (std::is_same_v<T, double>)
// [[nodiscard]] tf::vec3<T> cross_simd(const tf::vec3<T>& u, const tf::vec3<T>& v) noexcept {
//    const auto vec0 = _mm256_load_pd(&u[0]);
//    const auto vec1 = _mm256_load_pd(&v[0]);
//    const auto tmp0 = _mm256_shuffle_pd( vec0, vec0, _MM_SHUFFLE(3,0,2,1) );
//    const auto tmp1 = _mm256_shuffle_pd( vec1, vec1, _MM_SHUFFLE(3,1,0,2) );
//    const auto tmp2 = _mm256_mul_pd( tmp0, vec1 );
//    const auto tmp3 = _mm256_mul_pd( tmp0, tmp1 );
//    const auto tmp4 = _mm256_shuffle_pd( tmp2, tmp2, _MM_SHUFFLE(3,0,2,1) );
//    const auto tmp5 = _mm256_sub_pd( tmp3, tmp4 );
//    auto result =  tf::vec3<double>{};
//    _mm256_store_pd(&result[0], tmp5);
//    return result;
// }

// [[nodiscard]] inline tf::vec3<float> cross_simd_flt(const tf::vec3<float>& u, const tf::vec3<float>& v) noexcept {
//    const auto vec0 = _mm_load_ps(u.data);
//    const auto vec1 = _mm_load_ps(v.data);
//    const auto tmp0 = _mm_shuffle_ps( vec0, vec0, _MM_SHUFFLE(3,0,2,1) );
//    const auto tmp1 = _mm_shuffle_ps( vec1, vec1, _MM_SHUFFLE(3,1,0,2) );
//    const auto tmp2 = _mm_mul_ps( tmp0, vec1 );
//    const auto tmp3 = _mm_mul_ps( tmp0, tmp1 );
//    const auto tmp4 = _mm_shuffle_ps( tmp2, tmp2, _MM_SHUFFLE(3,0,2,1) );
//    return tf::vec3<float>{ _mm_sub_ps( tmp3, tmp4 ) };
// }

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
