//
// Created by cepheid on 10/17/24.
//

#ifndef ELECTROMAGNETICS_H
#define ELECTROMAGNETICS_H

#include "electromagnetics.param"
#include "em_definitions.h"
#include "bc_definitions.h"

namespace tf::electromagnetics
{
  //=================== Boundary Condition Selectors ===================
  //====================================================================
  template<typename T, EMFace F, EMSide S>
  using BCDataTL = TypeList<
    ReflectingData<T>,        // 0
    PeriodicData3D<T, F, S>,  // 1
    PmlData3D<T, F, S>        // 2
  >;

  template<typename T>
  using bcdata_t = boundaries::BCData<
    TypeListAt<SELECT_BCSOLVER[0], BCDataTL<T, EMFace::X, EMSide::Lo>>,
    TypeListAt<SELECT_BCSOLVER[1], BCDataTL<T, EMFace::X, EMSide::Hi>>,
    TypeListAt<SELECT_BCSOLVER[2], BCDataTL<T, EMFace::Y, EMSide::Lo>>,
    TypeListAt<SELECT_BCSOLVER[3], BCDataTL<T, EMFace::Y, EMSide::Hi>>,
    TypeListAt<SELECT_BCSOLVER[4], BCDataTL<T, EMFace::Z, EMSide::Lo>>,
    TypeListAt<SELECT_BCSOLVER[5], BCDataTL<T, EMFace::Z, EMSide::Hi>>
  >;

  template<EMFace F, EMSide S>
  using BCTypeTL = TypeList<
    ReflectingBC,     // 0
    Periodic3D<F, S>, // 1
    Pml3D<F, S>       // 2
  >;

  template<size_t I, EMFace F, EMSide S>
  using BCType = TypeListAt<I, BCTypeTL<F, S>>;

  template<typename TL>
  struct boundary_t {
    using Ex = boundaries::Boundary<TypeListAt<0, TL>>;
    using Ey = boundaries::Boundary<TypeListAt<1, TL>>;
    using Ez = boundaries::Boundary<TypeListAt<2, TL>>;
    using Hx = boundaries::Boundary<TypeListAt<3, TL>>;
    using Hy = boundaries::Boundary<TypeListAt<4, TL>>;
    using Hz = boundaries::Boundary<TypeListAt<5, TL>>;
  };

  using bcx0_t = boundary_t<BCType<SELECT_BCSOLVER[0], EMFace::X, EMSide::Lo>>;
  using bcx1_t = boundary_t<BCType<SELECT_BCSOLVER[1], EMFace::X, EMSide::Hi>>;
  using bcy0_t = boundary_t<BCType<SELECT_BCSOLVER[2], EMFace::Y, EMSide::Lo>>;
  using bcy1_t = boundary_t<BCType<SELECT_BCSOLVER[3], EMFace::Y, EMSide::Hi>>;
  using bcz0_t = boundary_t<BCType<SELECT_BCSOLVER[4], EMFace::Z, EMSide::Lo>>;
  using bcz1_t = boundary_t<BCType<SELECT_BCSOLVER[5], EMFace::Z, EMSide::Hi>>;

  //=================== Electromagnetics Selectors ===================
  //==================================================================
  template<typename T>
  using EMDataTL = TypeList<
    emdataNone<T>, // 0
    emdata3D<T>    // 1
  >; // Typelist for choosing type of EMData

  template<typename T>
  using emdata_t = TypeListAt<SELECT_EMSOLVER, EMDataTL<T>>;

  template<typename T>
  using EMTypeTL = TypeList<
    EMNull<T>, // 0
    EM3D<T>    // 1
  >; // Typelist (of typelists) for choosing type of EM Solver

  template<typename T>
  using EMType = TypeListAt<SELECT_EMSOLVER, EMTypeTL<emdata_t<T>>>; // Selects desired typelist of integrators

  template<typename T> using EXI = TypeListAt<0, EMType<T>>;
  template<typename T> using EYI = TypeListAt<1, EMType<T>>;
  template<typename T> using EZI = TypeListAt<2, EMType<T>>;
  template<typename T> using HXI = TypeListAt<3, EMType<T>>;
  template<typename T> using HYI = TypeListAt<4, EMType<T>>;
  template<typename T> using HZI = TypeListAt<5, EMType<T>>;

  template<typename T>
  using EMSolver = Electromagnetics<
    EXI<T>, EYI<T>, EZI<T>,
    HXI<T>, HYI<T>, HZI<T>,
    bcx0_t, bcx1_t,
    bcy0_t, bcy1_t,
    bcz0_t, bcz1_t
  >; // Auto-selects full BC's per face and EMSolver per field component based on chosen settings.

} // end namespace tf::old_em
#endif //ELECTROMAGNETICS_H
