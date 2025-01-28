#include <x86intrin.h>
#include <iostream>
#include <chrono>
#include <cassert>
#include <cmath>
#include <array>

// This function checks if two values are within tolerance of each other
template<class T>
bool approximatelyEqual(T a, T b, T tol = static_cast<T>(1e-15)) { return std::abs(a - b) <= tol; }

struct alignas(32) vec3 {
  double e[4]{};
};

std::ostream& operator<<(std::ostream& os, const vec3& v) {
  return os << "(" << v.e[0] << ", " << v.e[1] << ", " << v.e[2] << ", " << v.e[3] << ")";
}

// Method 5: Fewer swizzles, swizzle before subtraction
[[nodiscard]] inline static __m256d cross_product5(__m256d const& vec0, __m256d const& vec1) noexcept {
  const __m256d tmp0 = _mm256_permute4x64_pd(vec0, _MM_SHUFFLE(3,0,2,1) );
  const __m256d tmp1 = _mm256_permute4x64_pd(vec1, _MM_SHUFFLE(3,1,0,2) );
  const __m256d tmp2 = _mm256_mul_pd(tmp0, vec1 );
  const __m256d tmp3 = _mm256_mul_pd(tmp0, tmp1 );
  const __m256d tmp4 = _mm256_permute4x64_pd(tmp2, _MM_SHUFFLE(3,0,2,1) );
  return _mm256_sub_pd( tmp3, tmp4 );
}

vec3 cross(const vec3& u, const vec3& v) {
  // Performs u x v
  return {u.e[1] * v.e[2] - u.e[2] * v.e[1],
          u.e[2] * v.e[0] - u.e[0] * v.e[2],
          u.e[0] * v.e[1] - u.e[1] * v.e[0],
          0.0};
}

#define NUM 10000

void validate(const std::array<vec3, NUM>& test, const std::array<vec3, NUM>& solution) {
  for (std::size_t i = 0; i < NUM - 1; ++i) {
  //   std::cout << test[i] << " | " << solution[i] << std::endl;
  assert(approximatelyEqual(test[i].e[0], solution[i].e[0]));
  assert(approximatelyEqual(test[i].e[1], solution[i].e[1]));
  assert(approximatelyEqual(test[i].e[2], solution[i].e[2]));
  assert(approximatelyEqual(test[i].e[3], solution[i].e[3]));
  }
}

int main() {
  std::array<vec3, NUM> array{};
  std::array<vec3, NUM> i1_result{};
  std::array<vec3, NUM> reg_result{};

  for (size_t i = 0; i < NUM; ++i) {
    const double r1 = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
    const double r2 = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
    const double r3 = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
    array[i].e[0] = r1;
    array[i].e[1] = r2;
    array[i].e[2] = r3;
    array[i].e[3] = 0.0f;
  }

  const auto reg_start = std::chrono::high_resolution_clock::now();
  for (std::size_t i = 0; i < NUM - 1; ++i) {
    reg_result[i] = cross(array[i], array[i + 1]);
  }
  const auto reg_stop = std::chrono::high_resolution_clock::now() - reg_start;

  const auto i5_start = std::chrono::high_resolution_clock::now();
  for (std::size_t i = 0; i < NUM - 1; ++i) {
    const __m256d vec0 = _mm256_load_pd(array[i].e);
    const __m256d vec1 = _mm256_load_pd(array[i + 1].e);
    const __m256d val = cross_product5(vec0, vec1);
    _mm256_store_pd(i1_result[i].e, val);
  }
  const auto i5_stop = std::chrono::high_resolution_clock::now() - i5_start;
  validate(i1_result, reg_result);

  const auto i5_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(i5_stop);
  const auto reg_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(reg_stop);

  const auto i5df = static_cast<double>(i5_duration.count());
  const auto regf = static_cast<double>(reg_duration.count());

  const auto i5ave = static_cast<double>(i5df) / static_cast<double>(NUM - 1);
  const auto rave = static_cast<double>(regf) / static_cast<double>(NUM - 1);
  // const auto speedup = iave / rave * 100.0;

  std::cout << "Intrinics 5: " << i5_duration.count() << " ns (" << i5ave << " ns/op)" << std::endl;
  std::cout << "    Regular: " << reg_duration.count() << " ns (" << rave << " ns/op)" << std::endl;

  return 0;
}
