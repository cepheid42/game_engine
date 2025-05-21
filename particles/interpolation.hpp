#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

#include "program_params.hpp"
#include "math_utils.hpp"
#include "dbg.h"

#include <cmath>
#include <type_traits>

namespace tf::interp
{
template<int D>
constexpr auto rotateOrigin(const auto& p)
{
   using ret_type = std::remove_cvref_t<decltype(p)>;
   if constexpr (D == 0)
   {
      return ret_type{p[1], p[2], p[0]};
   }
   else if constexpr (D == 1)
   {
      return ret_type{p[2], p[0], p[1]};
   }
   else
   {
      return p;
   }
}

struct TSC
{
   static constexpr int         begin   = -1;
   static constexpr int         end     = 1;
   static constexpr std::size_t support = 3;

   static constexpr auto innerRadius(const compute_t x)
   {
      const auto val = 0.75_fp - math::SQR(x);
      return val; // < 1.0e-18 ? 0.0 : val;
   }

   static constexpr auto outerRadius(const compute_t x)
   {
      const auto val = 0.5_fp * math::SQR(1.5_fp - x);
      return val; // < 1.0e-18 ? 0.0 : val;
   }

   static constexpr auto operator()(const compute_t x)
   {
      const auto absx = std::abs(x);
      if (absx < 0.5_fp)
      {
         return innerRadius(x);
      }
      return outerRadius(absx);
   }
};

struct CIC
{
   static constexpr int         begin   = 0;
   static constexpr int         end     = 1;
   static constexpr std::size_t support = 2;

   static constexpr auto operator()(const compute_t x)
   {
      return 1.0_fp - std::abs(x);
   }
};

template<typename ParticleAssignFunctor>
struct Jit
{
   const compute_t particle_position;

   constexpr explicit Jit(const compute_t pos)
   : particle_position(pos)
   {}

   constexpr auto operator()(const compute_t grid_point) const
   {
      return ParticleAssignFunctor()(grid_point - particle_position);
   }
};

template<typename ParticleAssignFunctor>
struct Cached
{
   const std::array<compute_t, 3> shapeArray;

   constexpr explicit Cached(const compute_t pos)
   : shapeArray(std::move(ParticleAssignFunctor().shapeArray(pos)))
   {}

   constexpr auto operator()(const int grid_point) const
   {
      return shapeArray[grid_point - ParticleAssignFunctor::begin];
   }
};
} // end namespace tf::interp

#endif //INTERPOLATION_HPP
