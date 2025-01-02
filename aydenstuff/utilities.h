//
// Created by cepheid on 9/14/23.
//

#ifndef TRIFORCE_MATH_UTILITIES_H
#define TRIFORCE_MATH_UTILITIES_H

#include <limits>
#include <cmath>
#include <random>
#include <algorithm>
#include <array>
#include <source_location>
#include <variant>

#include "vector.h"

namespace tf::utilities::math {
  template<typename T> constexpr T SQR(T x) { return x * x; }
  template<typename T> constexpr T CUBE(T x) { return x * x * x; }

  /*
   * a - first input value
   * b - second input value
   * epsilon - maximum difference allowed, e.g. std::numeric_limits<T>::epsilon()
   * abs_th - minimum finite value representable, e.g. std::numeric_limits<T>::min()
   */
  template<typename T>
  bool nearly_equal(T a, T b,
                    T epsilon = T(128) * std::numeric_limits<T>::epsilon(),
                    T abs_th = std::numeric_limits<T>::min())
  {
    // Brought to you by StackOverflow
    // https://stackoverflow.com/questions/4915462/how-should-i-do-floating-point-comparison
    assert(std::numeric_limits<T>::epsilon() <= epsilon);
    assert(epsilon < T(1.0));

    if (a == b) { return true; }
    auto diff = std::abs(a - b);
    auto norm = std::min((std::abs(a) + std::abs(b)), std::numeric_limits<T>::max());
//  auto norm = std::min(std::abs(a + b), std::numeric_limits<T>::max());

    return diff < std::max(abs_th, epsilon * norm);
  }

  // This function checks if two values are within tolerance of each other
  template<class T>
  bool approximatelyEqual(T a, T b, T tol = static_cast<T>(1e-10)) { return std::abs(a - b) <= tol; }

  // This function does the same for points
  template <std::floating_point fp>
  bool approximatelyEqual(vec2<fp> a, vec2<fp> b, fp tol = static_cast<fp>(1.0e-10)) {
    return ((std::abs(a[0] - b[0]) <= tol) && (std::abs(a[1] - b[1]) <= tol));
  }


  // This function create a uniform spread of points over a given range
  // This should mimic Numpy's linspace function exactly.
  template<typename T>
  std::vector<T> linspace(T start, T stop, size_t n_points, const bool endpoint=true) {
    std::vector<T> result(n_points);
    if (endpoint) {
      n_points -= 1;
      result[result.size() - 1] = stop;
    }
    auto delta = (stop - start) / static_cast<T>(n_points);
    T val = start;
    for (size_t i = 0; i < n_points; ++i) {
      result[i] = val;
      val += delta;
    }
    return result;
  }


  // This function creates a set of points that are exactly halfway between the points of the provided vector
  template <std::floating_point fp>
  std::vector<fp> meanspace(const std::vector<fp>& input)
  {
    // Initialize the output vector
    std::vector<fp> output(input.size()-1);
    // Iterate through, taking the average of each set of points to find the halfway point
    for (size_t i = 0; i < output.size(); i++) {
      output[i] = 0.5 * (input[i + 1] + input[i]);
    }
    // Return the staggered vector
    return output;
  }

  template <std::floating_point fp>
  std::vector<fp> diffspace(const std::vector<fp>& input) {
    // Initialize the output vector
    std::vector<fp> output(input.size()-1);
    // Iterate through, taking the difference of each set of points
    for (size_t i = 0; i < output.size(); i++) {
      output[i] = input[i+1] - input[i];
    }
    // Return the staggered vector
    return output;
  }

  // Function returns vector that is logarithmically spaced similar to np.logspace
  // https://stackoverflow.com/questions/21429294/is-there-something-like-numpy-logspace-in-c
  template <std::floating_point fp>
  std::vector<fp> logspace(const fp first, const fp last, const size_t n, const fp base = 10.0)
  {
    std::vector<fp> vals(n);

    double current_value = first;
    double step = (last - first) / static_cast<fp>(n - 1);

    for (size_t i = 0; i < n; i++) {
     vals[i] = std::pow(base, current_value);
     current_value += step;
    }
    return vals;
  } // end logspace()


  // RNG class definition
  template<std::floating_point fp, template<typename> class Distribution, typename... Args>
  class RandomNumberGenerator {
  public:
    explicit RandomNumberGenerator(std::decay_t<Args>... args) : distribution{args...}, generator(init_mt_64()){}

    auto operator()() {
      return distribution(generator);
    }

  private:
    Distribution<fp> distribution;
    std::mt19937_64 generator;

    static inline std::mt19937_64 init_mt_64() {
      std::array<int, 624> seed_data{};
      std::random_device r; // @akis: this may be problematic for portability later, be careful
      std::generate_n(seed_data.data(), seed_data.size(), std::ref(r));
      std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
      std::mt19937_64 gen(seq);
      return gen;
    }
  };

  template <std::floating_point fp>
  using uniform_rng = RandomNumberGenerator<fp, std::uniform_real_distribution, fp, fp>;
  
  template <std::floating_point fp>
  using normal_rng = RandomNumberGenerator<fp, std::normal_distribution, fp, fp>;

} // end namespace Math

namespace TOML {
  template<typename T>
  vec2<T> getVec2Values(const toml::array& arr) {
    assert(arr.size() == 2);

    auto x = arr[0].value<T>();
    auto y = arr[1].value<T>();
    if (!x || !y) {
      throw std::invalid_argument("Value missing when trying to load Vec2.");
    }
    return {x.value(), y.value()};
  }
  
  namespace detail {
    template <typename T>
    struct load_as_array : public std::false_type {};
    
    template <typename T>
    struct load_as_array<vec2<T>> : public std::true_type {};
    
    template <typename T>
    struct load_as_array<vec3<T>> : public std::true_type {};
  }
  
  template<typename T>
  std::vector<T> getArrayValues(const toml::array& arr)
  {
    std::vector<T> output;
    for (auto& opt : arr) {
      if constexpr (detail::load_as_array<T>()) {
        auto& vec = *opt.as_array();
        output.emplace_back( getVec2Values<typename T::value_t>(vec) );
      } else {
        std::optional<T> val = opt.value<T>(); // unwrap optional
        if (!val) {
          throw std::invalid_argument("When trying to load an array, a value was missing.");
        }
        output.push_back(val.value());
      }
    }
    return output;
  }

  // <AJK> 260824: Can this be generalized out?
  template<>
  std::vector<vec2<fptype>> getArrayValues<vec2<fptype>>(const toml::array& arr)
  {
    std::vector<vec2<fptype>> output;
    for (auto& opt : arr) {
      auto& vec = *opt.as_array();
      output.emplace_back( getVec2Values<fptype>(vec) );
    }
    return output;
  }
} // end namespace TOML


namespace algorithms {

  struct shuffleVector {
    std::random_device rd;
    std::mt19937 g;

    shuffleVector() : g(rd()) {}

    template<typename T>
    void operator() (std::vector<T> &vec, const size_t start, const size_t stop) {
      std::shuffle(vec.begin() + start, vec.begin() + stop, g);
    }
  };

} // end namespace algorithms



namespace Lambert {
  // todo: cleanup Lambert function
  // following code for computing principal branch of Lambert function from stack overflow
  // https://stackoverflow.com/questions/73760130/accurate-computation-of-principal-branch-of-the-lambert-w-function-with-standard

  float expf_scale_pos_normal(float a, int scale);

  float logf_pos_normal(float a);

  /*
    Compute the principal branch of the Lambert W function, W_0. The maximum
    error in the positive half-plane is 1.49874 ulps and the maximum error in
    the negative half-plane is 2.56002 ulps
  */
  float Lambert_w0f(float z) {
    const float em1_fact_0 = 0.625529587f; // exp(-1)_factor_0
    const float em1_fact_1 = 0.588108778f; // exp(-1)_factor_1
    const float qe1 = 2.71828183f / 4.0f; //  0x1.5bf0a8p-01 // exp(1)/4
    float e, w, num, den, rden, redz, y, r;

    if (std::isnan(z) || (z == INFINITY) || (z == 0.0f)) return z + z;
    if (fabsf(z) < 1.220703125e-4f) return fmaf(-z, z, z); // 0x1.0p-13
    redz = fmaf(em1_fact_0, em1_fact_1, z); // z + exp(-1)
    if (redz < 0.0625f) { // expansion at -(exp(-1))
      r = sqrtf(redz);
      w = -1.23046875f;  // -0x1.3b0000p+0
      w = fmaf(w, r, 2.17185670f); //  0x1.15ff66p+1
      w = fmaf(w, r, -2.19554094f); // -0x1.19077cp+1
      w = fmaf(w, r, 1.92107077f); //  0x1.ebcb4cp+0
      w = fmaf(w, r, -1.81141856f); // -0x1.cfb920p+0
      w = fmaf(w, r, 2.33162979f); //  0x1.2a72d8p+1
      w = fmaf(w, r, -1.00000000f); // -0x1.000000p+0
    } else {
      /* Compute initial approximation. Based on: Roberto Iacono and John
         Philip Boyd, "New approximations to the principal real-valued branch
         of the Lambert W function", Advances in Computational Mathematics,
         Vol. 43, No. 6, December 2017, pp. 1403-1436
      */
      y = fmaf(2.0f, sqrtf(fmaf(qe1, z, 0.25f)), 1.0f);
      y = logf_pos_normal(fmaf(1.15262585f, y, -0.15262585f) /
                          fmaf(0.45906518f, logf_pos_normal(y), 1.0f));
      w = fmaf(2.0390625f, y, -1.0f);

      /* perform Newton iterations to refine approximation to full accuracy */
      for (int i = 0; i < 3; i++) {
        e = expf_scale_pos_normal(w, -3); // 0.125f * expf (w);
        num = fmaf(w, e, -0.125f * z);
        den = fmaf(w, e, e);
        rden = 1.0f / den;
        w = fmaf(-num, rden, w);
      }
    }
    return w;
  }

  float uint32_as_float(uint32_t a) {
    float r;
    memcpy(&r, &a, sizeof r);
    return r;
  }

  uint32_t float_as_uint32(float a) {
    uint32_t r;
    memcpy(&r, &a, sizeof r);
    return r;
  }

  /* exp(a) * 2**scale; positive normal results only! Maximum error 0.86565 ulp */
  float expf_scale_pos_normal(float a, int scale) {
    const float flt_int_cvt = 12582912.0f; // 0x1.8p23
    float f, r, j, t;
    uint32_t i;

    /* exp(a) = 2**i * exp(f); i = rintf (a / log(2)) */
    j = fmaf(1.442695f, a, flt_int_cvt); // // 0x1.715476p0 // log2(e)
    t = j - flt_int_cvt;
    f = fmaf(t, -6.93145752e-1f, a); // -0x1.62e400p-1  // log_2_hi
    f = fmaf(t, -1.42860677e-6f, f); // -0x1.7f7d1cp-20 // log_2_lo
    i = float_as_uint32(j);

    /* approximate r = exp(f) on interval [-log(2)/2, +log(2)/2] */
    r = 1.37805939e-3f;  // 0x1.694000p-10
    r = fmaf(r, f, 8.37312452e-3f); // 0x1.125edcp-7
    r = fmaf(r, f, 4.16695364e-2f); // 0x1.555b5ap-5
    r = fmaf(r, f, 1.66664720e-1f); // 0x1.555450p-3
    r = fmaf(r, f, 4.99999851e-1f); // 0x1.fffff6p-2
    r = fmaf(r, f, 1.00000000e+0f); // 0x1.000000p+0
    r = fmaf(r, f, 1.00000000e+0f); // 0x1.000000p+0

    /* exp(a) = 2**(i+scale) * r; */
    r = uint32_as_float(((i + scale) << 23) + float_as_uint32(r));
    return r;
  }

  /* compute natural logarithm of positive normals; maximum error: 0.85089 ulp */
  float logf_pos_normal(float a) {
    const float ln2 = 0.693147182f; // 0x1.62e430p-1 // log(2)
    const float two_to_m23 = 1.19209290e-7f; // 0x1.0p-23
    float m, r, s, t, i, f;
    int32_t e;

    /* log(a) = log(m * 2**i) = i * log(2) + log(m) */
    e = (float_as_uint32(a) - float_as_uint32(0.666666667f)) & 0xff800000;
    m = uint32_as_float(float_as_uint32(a) - e);
    i = (float) e * two_to_m23;

    /* log(m) = log1p(f) */
    f = m - 1.0f;
    s = f * f;

    /* compute log1p(f) for f in [-1/3, 1/3] */
    r = -0.130310059f;  // -0x1.0ae000p-3
    t = 0.140869141f;  //  0x1.208000p-3
    r = fmaf(r, s, -0.121483363f); // -0x1.f1988ap-4
    t = fmaf(t, s, 0.139814854f); //  0x1.1e5740p-3
    r = fmaf(r, s, -0.166846141f); // -0x1.55b36ep-3
    t = fmaf(t, s, 0.200120345f); //  0x1.99d8b2p-3
    r = fmaf(r, s, -0.249996200f); // -0x1.fffe02p-3
    r = fmaf(t, f, r);
    r = fmaf(r, f, 0.333331972f); //  0x1.5554fap-2
    r = fmaf(r, f, -0.500000000f); // -0x1.000000p-1
    r = fmaf(r, s, f);

    /* log(a) = i * log(2) + log(m) */
    r = fmaf(i, ln2, r);
    return r;
  }
} // end namespace Lambert


#endif //TRIFORCE_MATH_UTILITIES_H
