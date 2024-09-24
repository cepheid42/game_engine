//
// Created by akis on 9/9/24.
//

//#include <concepts>

//namespace scratch
//{
//  namespace concepts
//  {
//    template <typename T>
//    concept EMData = requires
//    {
//      requires array_like<T::Ex>;
//    };
//
//    template <typename T>
//    concept ElectricIntegrator = requires
//    {
//      T::advanceE;
//    };
//
//    template <typename T>
//    concept MagneticIntegrator = requires
//    {
//      T::advanceB;
//    };
//
//    template <typename T>
//    concept Solver = requires
//    {
//      requires EMData<typename T::emdata_t>;
//      T::advance;
//    };
//  }
//
//  using Mode1D = TypeList<ON, OFF, OFF,
//                          OFF, ON, OFF,
//                          ON, OFF, OFF>;
//
//  using ModeTE = ;
//  using ModeTM = ;
//  using Mode3d = ;
//
//
//
//
//  template <std::floating_point fp, typename T>
//  struct EMDataBase {
//    using ex_t = TypeListAt<0, T>;
//    ...
//
//    ex_t Ex; ey_t Ey; ez_t Ez;
//    bx_t Bx; by_t By; bz_t Bz;
//    jx_t Jx; jy_t Jy; jz_t Jz;
//  };
//
//  template <std::floating_point fp, concepts::Solver S>
//  struct BasicElectricIntegrator {
//    using emdata_t = S::emdata_t;
//
//    void advanceE(void) { return; }
//  };
//
//  template <std::floating_point fp, concepts::Solver S>
//  struct BasicMagneticIntegrator {
//    using emdata_t = S::emdata_t;
//
//    void advanceB(void) { return; }
//  };
//
//  template <std::floating_point fp, typename T>
//  struct EMSolver : public BasicElectricIntegrator<fp, EMSolver<fp>>, public BasicMagneticIntegrator<fp, EMSolver<fp>>
//  {
//    using array_t = void;
//
//    struct emdata_t : public EMDataBase<T> {
//      using typename EMDataBase<T>::ex_t;
//
//      ex_t ex_coeffs;
//      ...
//    };
//  };
//
////  template <std::floating_point fp, TypeList T>
////  struct ConformalSolver : public ... {
////    struct emdata_t : EMDataBase<>, ConformalDataBase <> {
////
////    };
////  };
//}

#include <iostream>
#include <cassert>

#include "../core/typelist.h"
#include "../aydenstuff/array.h"
#include "em_solver.h"

int main()
{
  using array_t = tf::types::Array2D<double>;
  using empty_t = tf::types::EmptyArray2D<double>;

  using emdata_tl = TypeList<
        empty_t,  // Ex
        empty_t,  // Ey
        array_t,  // Ez
        array_t,  // Hx
        array_t,  // Hy
        empty_t,  // Hz
        empty_t,  // Jx
        empty_t,  // Jy
        empty_t   // Jz
      >;

  return 0;
}