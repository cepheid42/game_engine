//
// Created by cepheid on 10/17/24.
//

#ifndef ELECTROMAGNETICS_H
#define ELECTROMAGNETICS_H

#include "em_definitions.h"

inline constexpr size_t SELECT_EMDATA = 1;
inline constexpr size_t SELECT_EMSOLVER = 1;

inline constexpr size_t SELECT_BCDATA = 1;
inline constexpr size_t SELECT_BCS[12] = {0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0};

using fp_t = double;
constexpr fp_t DIM = 1.0;
constexpr fp_t cfl = 0.95 / std::sqrt(DIM);
// constexpr size_t dPML = 10u;
// constexpr size_t nHalo = 2u;

template<typename T>
using EMDataTL = TypeList<emdataNone<T>, emdata1D<T>, emdataTM<T>, emdataTE<T>>;

template<typename T>
using emdata_t = TypeListAt<SELECT_EMDATA, EMDataTL<T>>;

template<typename T>
using EMTypeTL = TypeList<EMNull<T>, EM1D<T>, EMTM<T>, EMTE<T>>; // Typelist of typelists

template<typename T>
using EMType = TypeListAt<SELECT_EMSOLVER, EMTypeTL<emdata_t<T>>>; // Selects desired typelist of integrators

template<typename T>
using BCDataTL = TypeList<bcdataNone<T>, bcdata1d<T>, bcdataTM<T>, bcdataTE<T>>;

template<typename T>
using bcdata_t = TypeListAt<SELECT_BCDATA, BCDataTL<T>>;

template<typename T, bool HI, bool Forward>
using BCTypeTL = TypeList<NullBC<T>, Periodic1D<T>>;

template<size_t I, typename T, bool HI, bool Forward>
using BCType = TypeListAt<I, BCTypeTL<bcdata_t<T>, HI, Forward>>;

template<typename T>
using LoBoundaries = TypeList<
  TypeListAt<0, BCType<SELECT_BCS[0], T, false, false>>, // Ex
  TypeListAt<1, BCType<SELECT_BCS[1], T, false, false>>, // Ey
  TypeListAt<2, BCType<SELECT_BCS[2], T, false, false>>, // Ez
  TypeListAt<3, BCType<SELECT_BCS[3], T, false, false>>, // Hx
  TypeListAt<4, BCType<SELECT_BCS[4], T, false, false>>, // Hy
  TypeListAt<5, BCType<SELECT_BCS[5], T, false, false>>  // Hz
>;

template<typename T>
using HiBoundaries = TypeList<
  TypeListAt<0, BCType<SELECT_BCS[6], T, false, false>>,  // Ex
  TypeListAt<1, BCType<SELECT_BCS[7], T, false, false>>,  // Ey
  TypeListAt<2, BCType<SELECT_BCS[8], T, false, false>>,  // Ez
  TypeListAt<3, BCType<SELECT_BCS[9], T, false, false>>,  // Hx
  TypeListAt<4, BCType<SELECT_BCS[10], T, false, false>>, // Hy
  TypeListAt<5, BCType<SELECT_BCS[11], T, false, false>>  // Hz
>;

template<typename T>
using EMSolver = Electromagnetics<
  TypeListAt<0, EMType<T>>, // Ex
  TypeListAt<1, EMType<T>>, // Ey
  TypeListAt<2, EMType<T>>, // Ez
  TypeListAt<3, EMType<T>>, // Hx
  TypeListAt<4, EMType<T>>, // Hy
  TypeListAt<5, EMType<T>>, // Hz,
  LoBoundaries<T>,
  HiBoundaries<T>
>;

#endif //ELECTROMAGNETICS_H
