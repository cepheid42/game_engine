//
// Created by cepheid on 10/17/24.
//

#ifndef ELECTROMAGNETICS_H
#define ELECTROMAGNETICS_H

#include "em_definitions.h"

inline constexpr size_t SELECT_EMDATA = 4; // todo: these can be combined into one value?
inline constexpr size_t SELECT_EMSOLVER = 4;

inline constexpr size_t SELECT_BCDATA = 4;
inline constexpr size_t SELECT_BCS[2][6] = {
  {0, 0, 0, 0, 0, 0}, // Lo -> {Ex, Ey, Ez, Hx, Hy, Hz}
  {0, 0, 0, 0, 0, 0}  // Hi -> {Ex, Ey, Ez, Hx, Hy, Hz}
};

/*
 * Periodic2D in X -> Lo {0, 0, 0, 0, 2, 0}
 * Periodic2D in Y -> Lo {0, 0, 0, 2, 0, 0}
 */


using fp_t = double;
constexpr size_t DIM = 3;
constexpr fp_t cfl = 0.95 / std::sqrt(static_cast<fp_t>(DIM));


template<typename T>
using EMDataTL = TypeList<
  emdataNone<T>, // 0
  emdata1D<T>,   // 1
  emdataTM<T>,   // 2
  emdataTE<T>,   // 3
  emdata3D<T>    // 4
>; // Typelist for choosing type of EMData

template<typename T>
using emdata_t = TypeListAt<SELECT_EMDATA, EMDataTL<T>>;

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
using BCDataTL = TypeList<
  bcdataNone<T>, // 0
  bcdata1d<T>,   // 1
  bcdataTM<T>,   // 2
  bcdataTE<T>,   // 3
  bcdata3D<T>    // 4
>; // Typelist for choosing type of BCData

template<typename T>
using bcdata_t = TypeListAt<SELECT_BCDATA, BCDataTL<T>>;

template<typename T, typename Func>
using BCIntegratorTL = TypeList<BCIntegrator1D<T, Func>, BCIntegrator2D<T, Func>, BCIntegrator3D<T, Func>>;

template<size_t I, typename T, typename Func>
using BCIntegratorType = TypeListAt<I, BCIntegratorTL<bcdata_t<T>, Func>>;

template<typename T, EMFace Face>
using BCTypeTL = TypeList<
  ReflectingBC<typename bcdata_t<T>::array_t>, // 0
  Periodic1D<bcdata_t<T>, Face>,               // 1
  Periodic2D<bcdata_t<T>, Face>,               // 2
  Periodic3D<bcdata_t<T>, Face>                // 3
>; // Typelist of different Boundary Conditions (per field per face)

template<typename T, typename B, size_t I>
using BoundaryType = BCApplicator<
  B,
  BCIntegratorType<DIM - 1, T, TypeListAt<SELECT_BCS[I][0], BCTypeTL<T, B::face>>>,
  BCIntegratorType<DIM - 1, T, TypeListAt<SELECT_BCS[I][1], BCTypeTL<T, B::face>>>,
  BCIntegratorType<DIM - 1, T, TypeListAt<SELECT_BCS[I][2], BCTypeTL<T, B::face>>>,
  BCIntegratorType<DIM - 1, T, TypeListAt<SELECT_BCS[I][3], BCTypeTL<T, B::face>>>,
  BCIntegratorType<DIM - 1, T, TypeListAt<SELECT_BCS[I][4], BCTypeTL<T, B::face>>>,
  BCIntegratorType<DIM - 1, T, TypeListAt<SELECT_BCS[I][5], BCTypeTL<T, B::face>>>
>; // Auto-selects Boundaries based on face and SELECT_BCS settings

template<typename T>
using EMSolver = Electromagnetics<
  TypeListAt<0, EMType<T>>, // Ex
  TypeListAt<1, EMType<T>>, // Ey
  TypeListAt<2, EMType<T>>, // Ez
  TypeListAt<3, EMType<T>>, // Hx
  TypeListAt<4, EMType<T>>, // Hy
  TypeListAt<5, EMType<T>>, // Hz,
  BoundaryType<T, XLo, 0>,
  BoundaryType<T, YLo, 0>
  // BoundaryType<T, ZLo, 2>
>; // Auto-selects full BC's per face and EMSolver per field component based on chosen settings.

#endif //ELECTROMAGNETICS_H
