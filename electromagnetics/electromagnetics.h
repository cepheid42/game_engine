//
// Created by cepheid on 10/17/24.
//

#ifndef ELECTROMAGNETICS_H
#define ELECTROMAGNETICS_H

#include "em_definitions.h"

inline constexpr size_t SELECT_EMDATA = 1;
inline constexpr size_t SELECT_EMSOLVER = 1;

inline constexpr size_t SELECT_BCDATA = 1;
inline constexpr size_t SELECT_BCS[2][6] = {
  {0, 0, 0, 0, 1, 0},
  {0, 0, 0, 0, 0, 0}
};


//   {
//   1, // X0
//   0, // Y0
//   0, // Z0
//   0, // X1
//   0, // Y1
//   0  // Z1
// };


using fp_t = double;
constexpr size_t DIM = 1;
constexpr fp_t cfl = 0.95 / std::sqrt(static_cast<fp_t>(DIM));
// constexpr size_t dPML = 10u;
// constexpr size_t nHalo = 2u;


template<typename T>
using EMDataTL = TypeList<emdataNone<T>, emdata1D<T>, emdataTM<T>, emdataTE<T>>;

template<typename T>
using emdata_t = TypeListAt<SELECT_EMDATA, EMDataTL<T>>;

template<typename T>
using EMTypeTL = TypeList<EMNull<T>, EM1D<T>, EMTM<T>, EMTE<T>, EM3D<T>>; // Typelist of typelists

template<typename T>
using EMType = TypeListAt<SELECT_EMSOLVER, EMTypeTL<emdata_t<T>>>; // Selects desired typelist of integrators

template<typename T>
using BCDataTL = TypeList<bcdataNone<T>, bcdata1d<T>, bcdataTM<T>, bcdataTE<T>, bcdata3D<T>>;

template<typename T>
using bcdata_t = TypeListAt<SELECT_BCDATA, BCDataTL<T>>;

template<typename T, typename Func>
using BCIntegratorTL = TypeList<BCIntegrator1D<T, Func>, BCIntegrator2D<T, Func>, BCIntegrator3D<T, Func>>;

template<size_t I, typename T, typename Func>
using BCIntegratorType = TypeListAt<I, BCIntegratorTL<bcdata_t<T>, Func>>;

template<typename T>
using BCTypeTL = TypeList<
  // Reflecting1D<bcdata_t<T>>,
  // Reflecting2D<bcdata_t<T>>,
  // Reflecting3D<bcdata_t<T>>,
  ReflectingBC<typename bcdata_t<T>::array_t>,
  Periodic1D<bcdata_t<T>>,
  Periodic2D<bcdata_t<T>>,
  Periodic3D<bcdata_t<T>>
>;

template<typename T, typename B, size_t I>
using BoundaryType = BCApplicator<
  B,
  BCIntegratorType<DIM - 1, T, TypeListAt<SELECT_BCS[I][0], BCTypeTL<T>>>,
  BCIntegratorType<DIM - 1, T, TypeListAt<SELECT_BCS[I][1], BCTypeTL<T>>>,
  BCIntegratorType<DIM - 1, T, TypeListAt<SELECT_BCS[I][2], BCTypeTL<T>>>,
  BCIntegratorType<DIM - 1, T, TypeListAt<SELECT_BCS[I][3], BCTypeTL<T>>>,
  BCIntegratorType<DIM - 1, T, TypeListAt<SELECT_BCS[I][4], BCTypeTL<T>>>,
  BCIntegratorType<DIM - 1, T, TypeListAt<SELECT_BCS[I][5], BCTypeTL<T>>>
>;

template<typename T>
using EMSolver = Electromagnetics<
  TypeListAt<0, EMType<T>>, // Ex
  TypeListAt<1, EMType<T>>, // Ey
  TypeListAt<2, EMType<T>>, // Ez
  TypeListAt<3, EMType<T>>, // Hx
  TypeListAt<4, EMType<T>>, // Hy
  TypeListAt<5, EMType<T>>, // Hz,
  BoundaryType<T, XLo, 0>
  // BoundaryType<T, YLo, 1>,
  // BoundaryType<T, ZLo, 2>
>;

#endif //ELECTROMAGNETICS_H
