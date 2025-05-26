#ifndef EMSOLVER_CUH
#define EMSOLVER_CUH

#include <cuda_runtime.h>

// ===== Vector Types =====
// ========================
namespace tf::gpu
{
// ===== vec3 =====
// ================
template<typename T>
struct vec3
{
   __device__ constexpr vec3() = default;

   __device__ constexpr vec3(T e0, T e1, T e2) : data{e0, e1, e2} {}

   __device__ constexpr T&       operator[](std::size_t i) { return data[i]; }
   __device__ constexpr const T& operator[](std::size_t i) const { return data[i]; }

   __device__ [[nodiscard]] auto length_squared() const
   {
      auto sum = T(0.0);
      for (std::size_t i = 0; i < 3; i++) { sum += data[i] * data[i]; }
      return sum;
   }

   __device__ [[nodiscard]] T length() const { return std::hypot(data[0], data[1], data[2]); } //std::sqrt(length_squared()); }

   // Unary Negation
   __device__ vec3 operator-() const { return {-data[0], -data[1], -data[2]}; }

   __device__ vec3& operator+=(const vec3& v)
   {
      data[0] += v[0];
      data[1] += v[1];
      data[2] += v[2];
      return *this;
   }

   __device__ vec3& operator-=(const vec3& v)
   {
      data[0] -= v[0];
      data[1] -= v[1];
      data[2] -= v[2];
      return *this;
   }

   __device__ vec3& operator*=(const T s)
   {
      data[0] *= s;
      data[1] *= s;
      data[2] *= s;
      return *this;
   }

   __device__ vec3& operator/=(const T s)
   {
      data[0] /= s;
      data[1] /= s;
      data[2] /= s;
      return *this;
   }

   __device__ vec3& operator/=(const vec3& v)
   {
      data[0] /= v[0];
      data[1] /= v[1];
      data[2] /= v[2];
      return *this;
   }

   template<typename U>
   __device__ constexpr vec3<U> as_type() const
   {
      return {
         static_cast<U>(data[0]),
         static_cast<U>(data[1]),
         static_cast<U>(data[2])
      };
   }

   __device__ bool operator==(const vec3& v) const { return (data[0] == v[0] && data[1] == v[1] && data[2] == v[2]); }
   __device__ bool operator!=(const vec3& v) const { return !(*this == v); }

   template<std::size_t I>
   __device__ constexpr auto& get() & { return data[I]; }

   template<std::size_t I>
   __device__ constexpr const auto& get() const & { return data[I]; }

   T data[3];
}; // end struct tf::gpu::vec3
}  // end namespace tf

// ===== vec3-scalar Operators =====
// =================================
template<typename T>
__device__ tf::gpu::vec3<T> operator*(T s, const tf::gpu::vec3<T>& u)
{
   return {s * u[0], s * u[1], s * u[2]};
}

template<typename T>
__device__ tf::gpu::vec3<T> operator*(const tf::gpu::vec3<T>& u, T s)
{
   return s * u;
}

template<typename T>
__device__ tf::gpu::vec3<T> operator/(const tf::gpu::vec3<T>& u, T s)
{
   return (T(1) / s) * u;
}

template<typename T>
__device__ tf::gpu::vec3<T> operator+(const tf::gpu::vec3<T>& u, const T& v)
{
   return {u[0] + v, u[1] + v, u[2] + v};
}

template<typename T>
__device__ tf::gpu::vec3<T> operator+(const T& v, const tf::gpu::vec3<T>& u)
{
   return u + v;
}

template<typename T>
__device__ tf::gpu::vec3<T> operator-(const tf::gpu::vec3<T>& u, const T& v)
{
   return {u[0] - v, u[1] - v, u[2] - v};
}

template<typename T>
__device__ tf::gpu::vec3<T> operator-(const T& v, const tf::gpu::vec3<T>& u)
{
   return {v - u[0], v - u[1], v - u[2]};
}

// ===== vec3-vec3 Operators =====
// =================================
template<typename T>
__device__ tf::gpu::vec3<T> operator+(const tf::gpu::vec3<T>& u, const tf::gpu::vec3<T>& v)
{
   return {u[0] + v[0], u[1] + v[1], u[2] + v[2]};
}

template<typename T>
__device__ tf::gpu::vec3<T> operator-(const tf::gpu::vec3<T>& u, const tf::gpu::vec3<T>& v)
{
   return {u[0] - v[0], u[1] - v[1], u[2] - v[2]};
}

template<typename T>
__device__ tf::gpu::vec3<T> operator*(const tf::gpu::vec3<T>& u, const tf::gpu::vec3<T>& v)
{
   return {u[0] * v[0], u[1] * v[1], u[2] * v[2]};
}

template<typename T>
__device__ tf::gpu::vec3<T> operator/(const tf::gpu::vec3<T>& u, const tf::gpu::vec3<T>& v)
{
   return {u[0] / v[0], u[1] / v[1], u[2] / v[2]};
}

template<typename T>
__device__ tf::gpu::vec3<T> unit_vector(const tf::gpu::vec3<T>& u)
{
   return u / u.length();
}

template<typename T>
__device__ T dot(const tf::gpu::vec3<T>& u, const tf::gpu::vec3<T>& v)
{
   // Performs u @ v
   return (u[0] * v[0]) + (u[1] * v[1]) + (u[2] * v[2]);
}

template<typename T>
__device__ tf::gpu::vec3<T> cross(const tf::gpu::vec3<T>& u, const tf::gpu::vec3<T>& v)
{
   // Performs u x v
   return {
      u[1] * v[2] - u[2] * v[1],
      u[2] * v[0] - u[0] * v[2],
      u[0] * v[1] - u[1] * v[0]
   };
}

#endif //EMSOLVER_CUH
