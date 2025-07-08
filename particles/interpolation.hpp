#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

#include "program_params.hpp"
#include "math_utils.hpp"
// #include "dbg.h"

#include <cmath>

namespace tf::interp {
template<int D>
constexpr auto rotateOrigin(const auto x, const auto y, const auto z) {
   if constexpr (D == 0) {
      return vec3{y, z, x};
   }
   else if constexpr (D == 1) {
      return vec3{z, x, y};
   }
   else {
      return vec3{x, y, z};
   }
}

template<int D>
constexpr auto rotateOrigin(const auto& p) {
   if constexpr (D == 0) {
      return vec3{p[1], p[2], p[0]};
   }
   else if constexpr (D == 1) {
      return vec3{p[2], p[0], p[1]};
   }
   else {
      return p;
   }
}

struct NGP {
   static constexpr int         Begin   = 0;
   static constexpr int         End     = 0;
   static constexpr std::size_t Support = 1;

   static constexpr auto eval(const compute_t) {
      return 1.0;
   }

   static constexpr auto shape_array(const compute_t) {
      return std::array{1.0};
   }
};

struct CIC {
   static constexpr int         Begin   = 0;
   static constexpr int         End     = 1;
   static constexpr std::size_t Support = 2;

   static constexpr auto eval(const compute_t x) {
      return 1.0_fp - std::abs(x);
   }

   static constexpr auto shape_array(const compute_t x) {
      return std::array{eval(x - Begin), eval(x - End)};
   }
};

struct TSC {
   static constexpr int         Begin   = -1;
   static constexpr int         End     = 1;
   static constexpr std::size_t Support = 3;

   static constexpr auto innerRadius(const compute_t x) {
      return 0.75_fp - math::SQR(x);
   }

   static constexpr auto outerRadius(const compute_t x) {
      return 0.5_fp * math::SQR(1.5_fp - x);
   }

   static constexpr auto eval(const compute_t x) {
      const auto absx = std::abs(x);
      return absx <= 0.5_fp ? innerRadius(x) : outerRadius(absx);
   }

   static constexpr auto shape_array(const compute_t x) {
      return std::array{eval(x - Begin), eval(x), eval(x - End)};
   }
};


struct PQS {
   static constexpr int         Begin   = -1;
   static constexpr int         End     = 2;
   static constexpr std::size_t Support = 4;

   static constexpr auto innerRadius(const compute_t x) {
      return (2.0 / 3.0) - math::SQR(x) + 0.5 * math::CUBE(x);
   }

   static constexpr auto outerRadius(const compute_t x) {
      return (1.0 / 6.0) * math::CUBE(2.0 - x);
   }

   static constexpr auto eval(const compute_t x) {
      const auto absx = std::abs(x);
      return absx < 1.0_fp ? innerRadius(x) : outerRadius(absx);
   }
};

template<int ShapeOrder>
struct InterpolationShape;

template<>
struct InterpolationShape<0> {
   static constexpr auto Begin = NGP::Begin;
   static constexpr auto End = NGP::End;
   static constexpr auto Support = NGP::Support;
   static constexpr auto Order = Support - 1;
   
   using type = NGP;
};

template<>
struct InterpolationShape<1> {
   static constexpr auto Begin = CIC::Begin;
   static constexpr auto End = CIC::End;
   static constexpr auto Support = CIC::Support;
   static constexpr auto Order = Support - 1;
   
   using type = CIC;
};

template<>
struct InterpolationShape<2> {
   static constexpr auto Begin = TSC::Begin;
   static constexpr auto End = TSC::End;
   static constexpr auto Support = TSC::Support;
   static constexpr auto Order = Support - 1;
   
   using type = TSC;
};

template<>
struct InterpolationShape<3> {
   static constexpr auto Begin = PQS::Begin;
   static constexpr auto End = PQS::End;
   static constexpr auto Support = PQS::Support;
   static constexpr auto Order = Support - 1;

   using type = PQS;
};

template<typename Outer, typename Middle, typename Inner>
struct InterpolationStrategy {
   using OuterShape = Outer;
   using MiddleShape = Middle;
   using InnerShape = Inner;
};

// template<typename ParticleAssignFunctor>
// struct Jit {
//    static constexpr auto Begin = ParticleAssignFunctor::Begin;
//    static constexpr auto End = ParticleAssignFunctor::End;
//    static constexpr auto Support = ParticleAssignFunctor::Support;
//    static constexpr auto Order = Support - 1;
//
//    const compute_t particle_position;
//
//    constexpr explicit Jit(const compute_t pos)
//    : particle_position(pos)
//    {}
//
//    constexpr auto operator()(const compute_t offset) const {
//       return ParticleAssignFunctor::eval(particle_position - offset);
//    }
// };

// template<typename ParticleAssignFunctor>
// struct Cached {
//    // static constexpr int Begin = ParticleAssignFunctor::Begin;
//    // static constexpr int End = ParticleAssignFunctor::End;
//    // static constexpr int Support = ParticleAssignFunctor::Support;
//
//    const std::array<compute_t, 3> shapeArray;
//
//    constexpr explicit Cached(const compute_t pos)
//    : shapeArray(std::move(ParticleAssignFunctor().shapeArray(pos)))
//    {}
//
//    constexpr auto operator()(const int grid_point) const {
//       return shapeArray[grid_point - ParticleAssignFunctor::Begin];
//    }
// };
} // end namespace tf::interp

#endif //INTERPOLATION_HPP
