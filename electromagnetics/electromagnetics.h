//
// Created by cepheid on 10/17/24.
//

#ifndef ELECTROMAGNETICS_H
#define ELECTROMAGNETICS_H

#include "em_definitions.h"

// #define SELECT_EMDATA 1
// #define SELECT_EMSOLVER 1

inline constexpr size_t SELECT_EMDATA = 1;
inline constexpr size_t SELECT_BCDATA = 1;
inline constexpr size_t SELECT_BCS[12] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
inline constexpr size_t SELECT_EMSOLVER = 1;


template<typename T>
using EMDataTL = TypeList<emdataNone<T>, emdata1D<T>, emdataTM<T>, emdataTE<T>>;

template<typename T>
using emdata_t = TypeListAt<SELECT_EMDATA, EMDataTL<T>>;

template<typename T>
using BCDataTL = TypeList<bcdataNone<T>, bcdata1d<T>, bcdataTM<T>, bcdataTE<T>>;

template<typename T>
using bcdata_t = TypeListAt<SELECT_BCDATA, BCDataTL<T>>;

template<typename T, bool HI, bool Forward, bool Subtract>
using BCTypeTL = TypeList<NoneBC<T>, PeriodicBC<T>, PmlBC<T, HI, Forward, Subtract>>;

template<size_t I, typename T, bool HI, bool Forward, bool Subtract>
using BCType = TypeListAt<I, BCTypeTL<bcdata_t<T>, HI, Forward, Subtract>>;

template<typename T>
using EMTypeTL = TypeList<EMNull<T>, EM1D<T>, EMTM<T>, EMTE<T>>; // Typelist of typelists

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
  // TypeListAt<0, BCType<SELECT_BCS[0], T, false, false, >>,
>;

#endif //ELECTROMAGNETICS_H
