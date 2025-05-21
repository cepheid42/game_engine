#ifndef MATHS_H
#define MATHS_H

#include <cmath>
#include <vector>

namespace tf::math
{
template<typename T>
constexpr T SQR(T x) { return x * x; }

template<typename T>
constexpr T CUBE(T x) { return x * x * x; }

// This function create a uniform spread of points over a given range
// This should mimic Numpy's linspace function exactly.
template<typename T>
std::vector<T> linspace(T start, T stop, std::size_t n_points, const bool endpoint = true)
{
   std::vector<T> result(n_points);
   if (endpoint)
   {
      n_points -= 1;
      result[result.size() - 1] = stop;
   }
   auto delta = (stop - start) / static_cast<T>(n_points);
   T    val   = start;
   for (size_t i = 0; i < n_points; ++i)
   {
      result[i] = val;
      val += delta;
   }
   return result;
}

template<typename T>
std::size_t findIndex(const T loc, const std::vector<T>& vec, const bool right = false)
{
   if (loc <= vec[0]) { return 0; }
   const auto last = vec.size() - 1;
   if (loc >= vec[last]) { return last; }

   // Branchless Binary Search to find index for a given location on the full grid
   const T* base = vec.data();
   int      n    = vec.size();

   while (n > 1)
   {
      const auto half = n / 2;
      base            = (base[half] < loc) ? &base[half] : base;
      n -= half;
   }
   return (*base < loc) + base - vec.data() - 1 + right;
}

//  /*
//   * a - first input value
//   * b - second input value
//   * epsilon - maximum difference allowed, e.g. std::numeric_limits<T>::epsilon()
//   * abs_th - minimum finite value representable, e.g. std::numeric_limits<T>::min()
//   */
//  template<std::floating_point T>
//  bool nearly_equal(T a, T b,
//                    T epsilon = T(128.) * std::numeric_limits<T>::epsilon(),
//                    T abs_th = std::numeric_limits<T>::min())
//  {
//    // Brought to you by StackOverflow
//    // https://stackoverflow.com/questions/4915462/how-should-i-do-floating-point-comparison
//    assert(std::numeric_limits<T>::epsilon() <= epsilon);
//    assert(epsilon < T(1.0));
//
//    if (a == b) { return true; }
//    auto diff = std::abs(a - b);
//    auto norm = std::min((std::abs(a) + std::abs(b)), std::numeric_limits<T>::max());
////  auto norm = std::min(std::abs(a + b), std::numeric_limits<T>::max());
//
//    return diff < std::max(abs_th, epsilon * norm);
//  }
//
//  // This function checks if two values are within tolerance of each other
//  template<class T>
//  bool approximatelyEqual(T a, T b, T tol = static_cast<T>(1.0e-10)) { return std::abs(a - b) <= tol; }
//
//  // This function create a uniform spread of points over a given range
//  // This should mimic Numpy's linspace function exactly.
//  template<std::floating_point T>
//  std::vector<T> linspace(T start, T stop, size_t n_points, const bool endpoint=true) {
//    std::vector<T> result(n_points);
//    if (endpoint) {
//      n_points -= 1;
//      result[result.size() - 1] = stop;
//    }
//    auto delta = (stop - start) / static_cast<T>(n_points);
//    T val = start;
//    for (size_t i = 0; i < n_points; ++i) {
//      result[i] = val;
//      val += delta;
//    }
//    return result;
//  }
//
//
//  // This function creates a set of points that are exactly halfway between the points of the provided vector
//  template<std::floating_point T>
//  std::vector<T> meanspace(const std::vector<T>& input)
//  {
//    // Initialize the output vector
//    std::vector<T> output(input.size() - 1);
//    // Iterate through, taking the average of each set of points to find the halfway point
//    for (size_t i = 0u; i < output.size(); i++) {
//      output[i] = 0.5 * (input[i + 1] + input[i]);
//    }
//    // Return the staggered vector
//    return output;
//  }
//
//  template<std::floating_point T>
//  std::vector<T> diffspace(const std::vector<T>& input) {
//    // Initialize the output vector
//    std::vector<T> output(input.size() - 1);
//    // Iterate through, taking the difference of each set of points
//    for (size_t i = 0u; i < output.size(); i++) {
//      output[i] = input[i + 1] - input[i];
//    }
//    // Return the staggered vector
//    return output;
//  }
//
//  // Function returns vector that is logarithmically spaced similar to np.logspace
//  // https://stackoverflow.com/questions/21429294/is-there-something-like-numpy-logspace-in-c
//  template<std::floating_point T>
//  std::vector<T> logspace(const T first, const T last, const size_t n, const T base = 10.0)
//  {
//    std::vector<T> vals(n);
//
//    double current_value = first;
//    double step = (last - first) / static_cast<T>(n - 1);
//
//    for (size_t i = 0; i < n; i++) {
//     vals[i] = std::pow(base, current_value);
//     current_value += step;
//    }
//    return vals;
//  } // end logspace()
}

#endif //MATHS_H
