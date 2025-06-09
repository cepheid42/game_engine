#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

#include "program_params.hpp"
#include "math_utils.hpp"
// #include "dbg.h"

#include <cmath>
#include <type_traits>

namespace tf::interp {
template<int D>
constexpr auto rotateOrigin(const auto& x, const auto& y, const auto& z) {
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
   using ret_type = std::remove_cvref_t<decltype(p)>;
   if constexpr (D == 0) {
      return ret_type{p[1], p[2], p[0]};
   }
   else if constexpr (D == 1) {
      return ret_type{p[2], p[0], p[1]};
   }
   else {
      return p;
   }
}

struct TSC {
   static constexpr int         begin   = -1;
   static constexpr int         end     = 1;
   static constexpr std::size_t support = 3;

   static constexpr auto innerRadius(const compute_t x) {
      return 0.75_fp - math::SQR(x);
   }

   static constexpr auto outerRadius(const compute_t x) {
      return 0.5_fp * math::SQR(1.5_fp - x);
   }

   static constexpr auto operator()(const compute_t x) {
      const auto absx = std::abs(x);
      if (absx <= 0.5_fp) {
         return innerRadius(x);
      }
      return outerRadius(absx);
   }

   // static constexpr auto shapeArray(const compute_t x) {
   //    return std::array{
   //       outerRadius(std::abs(-0.5_fp - x)),
   //       innerRadius(0.5 - x),
   //       outerRadius(std::abs(1.5_fp - x))
   //    };
   // }
};

struct CIC {
   static constexpr int         begin   = -1;
   static constexpr int         end     = 0;
   static constexpr std::size_t support = 2;

   static constexpr auto operator()(const compute_t x) {
      return 1.0_fp - std::abs(x);
   }

   // static constexpr auto shapeArray(const compute_t x) {
   //    return std::array{
   //       1.0_fp - std::abs(-0.5_fp - x),
   //       1.0_fp - std::abs(0.5_fp - x),
   //       0.0
   //    };
   // }
};

template<typename ParticleAssignFunctor>
struct Jit {
   // static constexpr int begin = ParticleAssignFunctor::begin;
   // static constexpr int end = ParticleAssignFunctor::end;
   // static constexpr int support = ParticleAssignFunctor::support;

   const compute_t particle_position;

   constexpr explicit Jit(const compute_t pos)
   : particle_position(pos)
   {}

   constexpr auto operator()(const compute_t offset) const
      requires std::same_as<ParticleAssignFunctor, CIC>
   {
      return ParticleAssignFunctor()(offset - particle_position);
   }

   constexpr auto operator()(const compute_t offset) const
   requires std::same_as<ParticleAssignFunctor, TSC>
   {
      return ParticleAssignFunctor()(offset - particle_position + 0.5);
   }
};

// template<typename ParticleAssignFunctor>
// struct Cached {
//    // static constexpr int begin = ParticleAssignFunctor::begin;
//    // static constexpr int end = ParticleAssignFunctor::end;
//    // static constexpr int support = ParticleAssignFunctor::support;
//
//    const std::array<compute_t, 3> shapeArray;
//
//    constexpr explicit Cached(const compute_t pos)
//    : shapeArray(std::move(ParticleAssignFunctor().shapeArray(pos)))
//    {}
//
//    constexpr auto operator()(const int grid_point) const {
//       return shapeArray[grid_point - ParticleAssignFunctor::begin];
//    }
// };
} // end namespace tf::interp

#endif //INTERPOLATION_HPP
