#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

#include "program_params.hpp"
#include "math_utils.hpp"
#include "vec3.hpp"

#include <cmath>
#include <array>

namespace tf::interp
{
template<int D>
constexpr auto rotateOrigin(const auto x, const auto y, const auto z) {
   if      constexpr (D == 0) { return std::array{y, z, x}; } // 1, 2, 0
   else if constexpr (D == 1) { return std::array{z, x, y}; } // 2, 0, 1
   else                       { return std::array{x, y, z}; } // 0, 1, 2
}

// template<int D>
// constexpr auto rotateOrigin(const auto& p) {
//    if      constexpr (D == 0) { return vec3{p[1], p[2], p[0]}; } // y, z, x
//    else if constexpr (D == 1) { return vec3{p[2], p[0], p[1]}; } // z, x, y
//    else                       { return p; }                      // x, y, z
// }

struct NGP {
   static constexpr int         Begin   = -1;
   static constexpr int         End     = -1;
   static constexpr std::size_t Order   = 0;
   static constexpr std::size_t Support = 1;

   static constexpr auto eval(const double) {
      return 1.0;
   }

   static constexpr auto shape_array(const double) {
      return std::array{1.0};
   }

   static constexpr auto ds_array(const auto, const auto&) {
      return std::array{1.0};
   }
};

struct CIC {
   static constexpr int         Begin   = -1;
   static constexpr int         End     = 0;
   static constexpr std::size_t Order   = 1;
   static constexpr std::size_t Support = 2;

   static constexpr auto eval(const double x) {
      return 1.0 - std::abs(x);
   }

   static constexpr auto shape_array(const double x) {
      return std::array{eval(x - Begin), eval(x)};
   }

   static constexpr auto ds_array(const auto x1, const auto& s0) {
      return std::array{eval(x1 - Begin) - s0[0],
                        eval(x1)         - s0[1]};
   }
};

struct TSC {
   static constexpr int         Begin   = -1;
   static constexpr int         End     = 1;
   static constexpr std::size_t Order   = 2;
   static constexpr std::size_t Support = 3;

   static constexpr auto innerRadius(const auto x) {
      return 0.75 - math::SQR(x);
   }

   static constexpr auto outerRadius(const auto x) {
      return 0.5 * math::SQR(1.5 - x);
   }

   static constexpr auto eval(const auto x) {
      const auto absx = std::abs(x);
      return absx <= 0.5 ? innerRadius(absx) : outerRadius(absx);
   }

   static constexpr auto shape_array(const auto x) {
      return std::array{eval(x - Begin), eval(x), eval(x - End)};
   }

   static constexpr auto ds_array(const auto x1, const auto& s0) {
      return std::array{eval(x1 - Begin) - s0[0],
                        eval(x1)         - s0[1],
                        eval(x1 - End)   - s0[2]};
   }
};

struct PQS {
   static constexpr int         Begin   = -1;
   static constexpr int         End     = 2;
   static constexpr std::size_t Order   = 3;
   static constexpr std::size_t Support = 4;

   static constexpr auto innerRadius(const double x) {
      return (2.0 / 3.0) - math::SQR(x) + 0.5 * math::CUBE(x);
   }

   static constexpr auto outerRadius(const double x) {
      return (1.0 / 6.0) * math::CUBE(2.0 - x);
   }

   static constexpr auto eval(const double x) {
      const auto absx = std::abs(x);
      return absx <= 1.0 ? innerRadius(x) : outerRadius(absx);
   }

   static constexpr auto shape_array(const double x) {
      return std::array{eval(x - Begin), eval(x), eval(x - 1.0), eval(x - End)};
   }

   static constexpr auto ds_array(const auto x1, const auto& s0) {
      return std::array{eval(x1 - Begin) - s0[0],
                        eval(x1)         - s0[1],
                        eval(x1 - 1.0)   - s0[2],
                        eval(x1 - End)   - s0[2]};
   }
};

template<int ShapeOrder> struct InterpolationShape;
template<> struct InterpolationShape<0> { using Type = NGP; };
template<> struct InterpolationShape<1> { using Type = CIC; };
template<> struct InterpolationShape<2> { using Type = TSC; };
template<> struct InterpolationShape<3> { using Type = PQS; };

template<typename Outer, typename Middle, typename Inner>
struct InterpolationStrategy {
   using OuterShape = Outer;
   using MiddleShape = Middle;
   using InnerShape = Inner;
};
} // end namespace tf::interp
#endif //INTERPOLATION_HPP
