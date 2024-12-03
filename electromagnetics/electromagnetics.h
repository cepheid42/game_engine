//
// Created by cepheid on 10/17/24.
//

#ifndef ELECTROMAGNETICS_H
#define ELECTROMAGNETICS_H

#include "em_definitions.h"
#include "bc_definitions.h"

inline constexpr size_t SELECT_EMSOLVER = 4;
inline constexpr size_t SELECT_BCSOLVER[6] = {8, 8, 8, 8, 8, 8}; // Xlo, Xhi, Ylo, Yhi, Zlo, Zhi

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
  TypeListAt<SELECT_BCSOLVER[0], BCDataTL<T, EMFace::X, EMSide::Lo>>,
  TypeListAt<SELECT_BCSOLVER[1], BCDataTL<T, EMFace::X, EMSide::Hi>>,
  TypeListAt<SELECT_BCSOLVER[2], BCDataTL<T, EMFace::Y, EMSide::Lo>>,
  TypeListAt<SELECT_BCSOLVER[3], BCDataTL<T, EMFace::Y, EMSide::Hi>>,
  TypeListAt<SELECT_BCSOLVER[4], BCDataTL<T, EMFace::Z, EMSide::Lo>>,
  TypeListAt<SELECT_BCSOLVER[5], BCDataTL<T, EMFace::Z, EMSide::Hi>>
>;

// template<EMFace F, EMSide S, bool Negate>
// using BCTypeTL = TypeList<
//   ReflectingBC,        // 0
//   Periodic1D,          // 1
//   Periodic2D<F, S>,    // 2
//   Periodic2D<F, S>,    // 3 (doubled for TM/TE so the numbers match)
//   Periodic3D<F, S>,    // 4
//   Pml1D<S>,            // 5
//   Pml2D<F, S, Negate>, // 6
//   Pml2D<F, S, Negate>, // 7 (doubled for TM/TE so the numbers match)
//   Pml3D<F, S, Negate>  // 8
// >;

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
using EMSolver = Electromagnetics<
  TypeListAt<0, EMType<T>>, // Ex
  TypeListAt<1, EMType<T>>, // Ey
  TypeListAt<2, EMType<T>>, // Ez
  TypeListAt<3, EMType<T>>, // Hx
  TypeListAt<4, EMType<T>>, // Hy
  TypeListAt<5, EMType<T>>  // Hz,
>; // Auto-selects full BC's per face and EMSolver per field component based on chosen settings.

#endif //ELECTROMAGNETICS_H
