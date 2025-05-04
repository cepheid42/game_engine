#ifndef TRIFORCE_VECTOR_H
#define TRIFORCE_VECTOR_H

#include <cmath>
#include <sstream>
#include <string>
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
    
    [[nodiscard]] T length() const { return std::hypot(data[0], data[1], data[2]); }//std::sqrt(length_squared()); }

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
    constexpr vec3<U> as_type() const {
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

// [[nodiscard]] inline tf::vec3<double> cross_simd_dub(const tf::vec3<double>& u, const tf::vec3<double>& v) noexcept {
//   const auto vec0 = _mm256_load_pd(u.data);
//   const auto vec1 = _mm256_load_pd(v.data);
//   const auto tmp0 = _mm256_shuffle_pd( vec0, vec0, _MM_SHUFFLE(3,0,2,1) );
//   const auto tmp1 = _mm256_shuffle_pd( vec1, vec1, _MM_SHUFFLE(3,1,0,2) );
//   const auto tmp2 = _mm256_mul_pd( tmp0, vec1 );
//   const auto tmp3 = _mm256_mul_pd( tmp0, tmp1 );
//   const auto tmp4 = _mm256_shuffle_pd( tmp2, tmp2, _MM_SHUFFLE(3,0,2,1) );
//   return tf::vec3<double>{ _mm256_sub_pd( tmp3, tmp4 ) };
// }
//
// [[nodiscard]] inline tf::vec3<float> cross_simd_flt(const tf::vec3<float>& u, const tf::vec3<float>& v) noexcept {
//   const auto vec0 = _mm_load_ps(u.data);
//   const auto vec1 = _mm_load_ps(v.data);
//   const auto tmp0 = _mm_shuffle_ps( vec0, vec0, _MM_SHUFFLE(3,0,2,1) );
//   const auto tmp1 = _mm_shuffle_ps( vec1, vec1, _MM_SHUFFLE(3,1,0,2) );
//   const auto tmp2 = _mm_mul_ps( tmp0, vec1 );
//   const auto tmp3 = _mm_mul_ps( tmp0, tmp1 );
//   const auto tmp4 = _mm_shuffle_ps( tmp2, tmp2, _MM_SHUFFLE(3,0,2,1) );
//   return tf::vec3<float>{ _mm_sub_ps( tmp3, tmp4 ) };
// }

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
