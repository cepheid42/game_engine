//
// Created by cepheid on 10/17/24.
//

#ifndef ELECTROMAGNETICS_H
#define ELECTROMAGNETICS_H

#include "em_definitions.h"

inline constexpr size_t SELECT_EMDATA = 2;
inline constexpr size_t SELECT_EMSOLVER = 2;

inline constexpr size_t SELECT_BCDATA = 2;
inline constexpr size_t SELECT_BCInt[6] = {1, 0, 0, 0, 0, 0};
inline constexpr size_t SELECT_BCS[6] = {
  1, // X0
  0, // Y0
  0, // Z0
  0, // X1
  0, // Y1
  0, // Z1
};


using fp_t = double;
constexpr fp_t DIM = 2.0;
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

template<typename T, bool HI, typename... IDXS>
using BCTypeTL = TypeList<NoneBC<T>, PeriodicBC<T, IDXS...>>;

template<size_t I, typename T, bool HI, typename... IDXS>
using BCType = TypeListAt<I, BCTypeTL<bcdata_t<T>, HI, IDXS...>>;

template<typename T, typename UpdateFunc>
using BCIntegratorTL = TypeList<BCIntegratorNull<T>, BCIntegrator1D<T, UpdateFunc>, BCIntegrator2D<T, UpdateFunc>>;

template<size_t I, typename T, typename UpdateFunc>
using BCIntegratorType = TypeListAt<I, BCIntegratorTL<bcdata_t<T>, UpdateFunc>>;

template<typename T, typename... IDXS>
using LoBoundaries = TypeList<
  BCIntegratorType<SELECT_BCInt[0], T, BCType<SELECT_BCS[0], T, false, IDXS...>>,
  BCIntegratorType<SELECT_BCInt[1], T, BCType<SELECT_BCS[0], T, false, IDXS...>>,
  BCIntegratorType<SELECT_BCInt[2], T, BCType<SELECT_BCS[0], T, false, IDXS...>>
>;

template<typename T>
using HiBoundaries = TypeList<
  // TypeListAt<0, BCType<SELECT_BCS[3], T, IDXS...>>,
  // TypeListAt<1, BCType<SELECT_BCS[4], T, IDXS...>>,
  // TypeListAt<2, BCType<SELECT_BCS[5], T, IDXS...>>
>;

template<typename T, typename... IDXS>
using EMSolver = Electromagnetics<
  TypeListAt<0, EMType<T>>, // Ex
  TypeListAt<1, EMType<T>>, // Ey
  TypeListAt<2, EMType<T>>, // Ez
  TypeListAt<3, EMType<T>>, // Hx
  TypeListAt<4, EMType<T>>, // Hy
  TypeListAt<5, EMType<T>>, // Hz,
  LoBoundaries<T, IDXS...>
  // HiBoundaries<T>
>;

#endif //ELECTROMAGNETICS_H
