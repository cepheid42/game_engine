//
// Created by cepheid on 10/17/24.
//

#ifndef ELECTROMAGNETICS_H
#define ELECTROMAGNETICS_H

#include "em_definitions.h"
#include "bc_definitions.h"

inline constexpr size_t SELECT_EMSOLVER = 4;
inline constexpr size_t SELECT_BCDATA[6] = {8, 8, 8, 8, 8, 8}; // Xlo, Xhi, Ylo, Yhi, Zlo, Zhi
inline constexpr size_t SELECT_BCSOLVER[6] = {9, 9, 10, 10, 11, 11}; // Xlo, Xhi, Ylo, Yhi, Zlo, Zhi

//=================== Boundary Condition Selectors ===================
//====================================================================
template<typename T, EMFace F, EMSide S>
using BCDataTL = TypeList<
  ReflectingFace<T, S>,     // 0
  PeriodicFace1D<T, F, S>,  // 1
  PeriodicFaceTM<T, F, S>,  // 2
  PeriodicFaceTE<T, F, S>,  // 3
  PeriodicFace3D<T, F, S>,  // 4
  PmlFace1D<T, F, S>,       // 5
  PmlFaceTM<T, F, S>,       // 6
  PmlFaceTE<T, F, S>,       // 7
  PmlFace3D<T, F, S>        // 8
>;

template<typename T>
using bcdata_t = BCData<
  TypeListAt<SELECT_BCDATA[0], BCDataTL<T, EMFace::X, EMSide::Lo>>,
  TypeListAt<SELECT_BCDATA[1], BCDataTL<T, EMFace::X, EMSide::Hi>>,
  TypeListAt<SELECT_BCDATA[2], BCDataTL<T, EMFace::Y, EMSide::Lo>>,
  TypeListAt<SELECT_BCDATA[3], BCDataTL<T, EMFace::Y, EMSide::Hi>>,
  TypeListAt<SELECT_BCDATA[4], BCDataTL<T, EMFace::Z, EMSide::Lo>>,
  TypeListAt<SELECT_BCDATA[5], BCDataTL<T, EMFace::Z, EMSide::Hi>>
>;

template<EMFace F, EMSide S>
using BCTypeTL = TypeList<
  BCNull,           // 0
  Periodic1D,       // 1
  Periodic2D<F, S>, // 2
  Periodic3D<F, S>, // 3
  PmlX1D<S>,        // 4
  PmlXTM<S>,        // 5
  PmlYTM<S>,        // 6
  PmlXTE<S>,        // 7
  PmlYTE<S>,        // 8
  PmlX3D<S>,        // 9
  PmlY3D<S>,        // 10
  PmlZ3D<S>         // 11
>;

template<size_t I, EMFace F, EMSide S>
using BCType = TypeListAt<I, BCTypeTL<F, S>>;

template<typename TL>
struct boundary_t {
  using Ex = Boundary<TypeListAt<0, TL>>;
  using Ey = Boundary<TypeListAt<1, TL>>;
  using Ez = Boundary<TypeListAt<2, TL>>;
  using Hx = Boundary<TypeListAt<3, TL>>;
  using Hy = Boundary<TypeListAt<4, TL>>;
  using Hz = Boundary<TypeListAt<5, TL>>;
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
  emdata1D<T>,   // 1
  emdataTM<T>,   // 2
  emdataTE<T>,   // 3
  emdata3D<T>    // 4
>; // Typelist for choosing type of EMData

template<typename T>
using emdata_t = TypeListAt<SELECT_EMSOLVER, EMDataTL<T>>;

template<typename T>
using EMTypeTL = TypeList<
  EMNull<T>, // 0
  EM1D<T>,   // 1
  EMTM<T>,   // 2
  EMTE<T>,   // 3
  EM3D<T>    // 4
>; // Typelist (of typelists) for choosing type of EM Solver

template<typename T>
using EMType = TypeListAt<SELECT_EMSOLVER, EMTypeTL<emdata_t<T>>>; // Selects desired typelist of integrators

template<typename T>
using EXI = TypeListAt<0, EMType<T>>;

template<typename T>
using EYI = TypeListAt<1, EMType<T>>;

template<typename T>
using EZI = TypeListAt<2, EMType<T>>;

template<typename T>
using HXI = TypeListAt<3, EMType<T>>;

template<typename T>
using HYI = TypeListAt<4, EMType<T>>;

template<typename T>
using HZI = TypeListAt<5, EMType<T>>;

template<typename T>
using EMSolver = Electromagnetics<
  EXI<T>, EYI<T>, EZI<T>,
  HXI<T>, HYI<T>, HZI<T>,
  bcx0_t, bcx1_t,
  bcy0_t, bcy1_t,
  bcz0_t, bcz1_t
>; // Auto-selects full BC's per face and EMSolver per field component based on chosen settings.

#endif //ELECTROMAGNETICS_H
