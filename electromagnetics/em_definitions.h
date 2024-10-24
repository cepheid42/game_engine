//
// Created by cepheid on 10/17/24.
//

#ifndef EM_DEFINITIONS_H
#define EM_DEFINITIONS_H

// #include "em_traits.h"
#include "em_data.h"
#include "em_solver.h"
#include "bc_data.h"
#include "boundaries.h"

//=================== EMData Definitions ========================
//===============================================================
template<typename Array>
struct enabled {
  using value_t = typename Array::value_t;
  using dimension_t = typename Array::dimension_t;
  using array_t = Array;
};

template<typename Array>
struct disabled {
  using value_t = typename Array::value_t;
  using dimension_t = typename Array::dimension_t;
  using array_t = Array;
};

template<typename T>
using emdataNone = EMData<
  disabled<EmptyArray1D<T>>, // Ex
  disabled<EmptyArray1D<T>>, // Ey
  disabled<EmptyArray1D<T>>, // Ez
  disabled<EmptyArray1D<T>>, // Hx
  disabled<EmptyArray1D<T>>, // Hy
  disabled<EmptyArray1D<T>>  // Hz
>;

template<typename T>
using emdata1D = EMData<
  disabled<EmptyArray1D<T>>, // Ex
  disabled<EmptyArray1D<T>>, // Ey
  enabled<Array1D<T>>,       // Ez
  disabled<EmptyArray1D<T>>, // Hx
  enabled<Array1D<T>>,       // Hy
  disabled<EmptyArray1D<T>>  // Hz
>;

template<typename T>
using emdataTM = EMData<
  disabled<EmptyArray2D<T>>, // Ex
  disabled<EmptyArray2D<T>>, // Ey
  enabled<Array2D<T>>,       // Ez
  enabled<Array2D<T>>,       // Hx
  enabled<Array2D<T>>,       // Hy
  disabled<EmptyArray2D<T>>  // Hz
>;

template<typename T>
using emdataTE = EMData<
  enabled<Array2D<T>>,       // Ex
  enabled<Array2D<T>>,       // Ey
  disabled<EmptyArray2D<T>>, // Ez
  disabled<EmptyArray2D<T>>, // Hx
  disabled<EmptyArray2D<T>>, // Hy
  enabled<Array2D<T>>        // Hz
>;

//=================== Electromagnetics Definitions ========================
//=========================================================================
template<typename T>
using EMNull = TypeList<
  /* Ex */ FieldIntegratorNull<T>,
  /* Ey */ FieldIntegratorNull<T>,
  /* Ez */ FieldIntegratorNull<T>,
  /* Hx */ FieldIntegratorNull<T>,
  /* Hy */ FieldIntegratorNull<T>,
  /* Hz */ FieldIntegratorNull<T>
>;

template<typename T>
using EM1D = TypeList<
  /* Ex */ FieldIntegratorNull<T>,
  /* Ey */ FieldIntegratorNull<T>,
  /* Ez */ FieldIntegrator1D<T, FieldUpdate<Derivative::DX, Derivative::NoOp, false, size_t>>,
  /* Hx */ FieldIntegratorNull<T>,
  /* Hy */ FieldIntegrator1D<T, FieldUpdate<Derivative::DX, Derivative::NoOp, true, size_t>>,
  /* Hz */ FieldIntegratorNull<T>
>;

template<typename T>
using EMTM = TypeList<
  /* Ex */ FieldIntegratorNull<T>,
  /* Ey */ FieldIntegratorNull<T>,
  /* Ez */ FieldIntegrator2D<T, FieldUpdate<Derivative::DX, Derivative::DY, false, size_t, size_t>>,
  /* Hx */ FieldIntegrator2D<T, FieldUpdate<Derivative::NoOp, Derivative::DY, true, size_t, size_t>>,
  /* Hy */ FieldIntegrator2D<T, FieldUpdate<Derivative::DX, Derivative::NoOp, true, size_t, size_t>>,
  /* Hz */ FieldIntegratorNull<T>
>;

template<typename T>
using EMTE = TypeList<
  /* Ex */ FieldIntegrator2D<T, FieldUpdate<Derivative::DY, Derivative::NoOp, false, size_t, size_t>>,
  /* Ey */ FieldIntegrator2D<T, FieldUpdate<Derivative::NoOp, Derivative::DX, false, size_t, size_t>>,
  /* Ez */ FieldIntegratorNull<T>,
  /* Hx */ FieldIntegratorNull<T>,
  /* Hy */ FieldIntegratorNull<T>,
  /* Hz */ FieldIntegrator2D<T, FieldUpdate<Derivative::DY, Derivative::DX, true, size_t, size_t>>
>;

//=================== BCData Definitions ========================
//===============================================================
template<typename T>
using bcdataNone = BCData<
  disabled<EmptyArray1D<T>>, // Ex
  disabled<EmptyArray1D<T>>, // Ey
  disabled<EmptyArray1D<T>>, // Ez
  disabled<EmptyArray1D<T>>, // Hx
  disabled<EmptyArray1D<T>>, // Hy
  disabled<EmptyArray1D<T>>  // Hz
>;

template<typename T>
using bcdata1d = BCData<
  disabled<EmptyArray1D<T>>, // Ex
  disabled<EmptyArray1D<T>>, // Ey
  enabled<Array1D<T>>,       // Ez
  disabled<EmptyArray1D<T>>, // Hx
  enabled<Array1D<T>>,       // Hy
  disabled<EmptyArray1D<T>>  // Hz
>;

template<typename T>
using bcdataTM = BCData<
  disabled<EmptyArray2D<T>>, // Ex
  disabled<EmptyArray2D<T>>, // Ey
  enabled<Array2D<T>>,       // Ez
  enabled<Array2D<T>>,       // Hx
  enabled<Array2D<T>>,       // Hy
  disabled<EmptyArray2D<T>>  // Hz
>;

template<typename T>
using bcdataTE = BCData<
  enabled<Array2D<T>>,       // Ex
  enabled<Array2D<T>>,       // Ey
  disabled<EmptyArray2D<T>>, // Ez
  disabled<EmptyArray2D<T>>, // Hx
  disabled<EmptyArray2D<T>>, // Hy
  enabled<Array2D<T>>        // Hz
>;

//=================== Boundary Conditions Definitions ========================
//=========================================================================
template<typename T>
using NullBC = TypeList<
  /* Ex */ BCIntegratorNull<T>,
  /* Ey */ BCIntegratorNull<T>,
  /* Ez */ BCIntegratorNull<T>,
  /* Hx */ BCIntegratorNull<T>,
  /* Hy */ BCIntegratorNull<T>,
  /* Hz */ BCIntegratorNull<T>
>;

template<typename T>
using Periodic1D = TypeList<
  /* Ex */ BCIntegratorNull<T>,
  /* Ey */ BCIntegratorNull<T>,
  /* Ez */ BCIntegrator1D<T, PeriodicBC<Array1D<T>, size_t>>,
  /* Hx */ BCIntegratorNull<T>,
  /* Hy */ BCIntegrator1D<T, PeriodicBC<Array1D<T>, size_t>>,
  /* Hz */ BCIntegratorNull<T>
>;

// template<typename T>
// using PeriodicTM = TypeList<
//   /* Ex */ BCIntegratorNull<T>,
//   /* Ey */ BCIntegratorNull<T>,
//   /* Ez */ BCIntegrator2D<T, PeriodicBC<Array2D<T>>>,
//   /* Hx */ BCIntegrator2D<T, PeriodicBC<Array2D<T>>>,
//   /* Hy */ BCIntegrator2D<T, PeriodicBC<Array2D<T>>>,
//   /* Hz */ BCIntegratorNull<T>
// >;

// template<typename T>
// using PML1D = TypeList<
//   /* Ex */ BCIntegratorNull<T>,
//   /* Ey */ BCIntegratorNull<T>,
//   /* Ez */ PmlBC<Array1D<T>>,
//   /* Hx */ BCIntegratorNull<T>,
//   /* Hy */ PeriodicBC<EmptyArray1D<T>>,
//   /* Hz */ BCIntegratorNull<T>
// >;

#endif //EM_DEFINITIONS_H
