//
// Created by cepheid on 10/17/24.
//

#ifndef EM_DEFINITIONS_H
#define EM_DEFINITIONS_H

#include "em_traits.h"
#include "em_data.h"
#include "em_solver.h"

//=================== EMData Definitions ========================
//===============================================================
template<typename Array>
struct enable_field {
  using value_t = typename Array::value_t;
  using dimension_t = typename Array::dimension_t;
  using array_t = Array;
};

template<typename Array>
struct disable_field {
  using value_t = typename Array::value_t;
  using dimension_t = typename Array::dimension_t;
  using array_t = EmptyArray<value_t, dimension_t::value>;
};

template<typename T>
using emdataNone = EMData<
  disable_field<EmptyArray1D<T>>, // Ex
  disable_field<EmptyArray1D<T>>, // Ey
  disable_field<EmptyArray1D<T>>, // Ez
  disable_field<EmptyArray1D<T>>, // Hx
  disable_field<EmptyArray1D<T>>, // Hy
  disable_field<EmptyArray1D<T>>  // Hz
>;


/*
 *
 * todo: Test changing integrator template to functor works correctly
 *
 */

template<typename T>
using emdata1D = EMData<
  disable_field<EmptyArray1D<T>>, // Ex
  disable_field<EmptyArray1D<T>>, // Ey
  enable_field<Array1D<T>>,       // Ez
  disable_field<Array1D<T>>,      // Hx
  enable_field<Array1D<T>>,       // Hy
  disable_field<Array1D<T>>       // Hz
>;

template<typename T>
using emdataTM = EMData<
  disable_field<EmptyArray2D<T>>, // Ex
  disable_field<EmptyArray2D<T>>, // Ey
  enable_field<Array2D<T>>,       // Ez
  enable_field<Array2D<T>>,       // Hx
  enable_field<Array2D<T>>,       // Hy
  disable_field<Array2D<T>>       // Hz
>;

template<typename T>
using emdataTE = EMData<
  enable_field<Array2D<T>>,       // Ex
  enable_field<Array2D<T>>,       // Ey
  disable_field<EmptyArray2D<T>>, // Ez
  disable_field<EmptyArray2D<T>>, // Hx
  disable_field<EmptyArray2D<T>>, // Hy
  enable_field<Array2D<T>>        // Hz
>;

//=================== Electromagnetics Definitions ========================
//=========================================================================
template<typename T>
// requires is_1D_fields<T, EmptyArray<typename T::value_t, T::dimension_t::value>>
using EMNull = TypeList<
  /* Ex */ FieldIntegratorNull<T>,
  /* Ey */ FieldIntegratorNull<T>,
  /* Ez */ FieldIntegratorNull<T>,
  /* Hx */ FieldIntegratorNull<T>,
  /* Hy */ FieldIntegratorNull<T>,
  /* Hz */ FieldIntegratorNull<T>
>;

template<typename T>
// requires is_1D_fields<T, EmptyArray<typename T::value_t, T::dimension_t::value>>
using EM1D = TypeList<
  /* Ex */ FieldIntegratorNull<T>,
  /* Ey */ FieldIntegratorNull<T>,
  /* Ez */ FieldIntegrator1D<T, FieldUpdate<Derivative::DX, Derivative::NoOp, false, size_t>>,
  /* Hx */ FieldIntegratorNull<T>,
  /* Hy */ FieldIntegrator1D<T, FieldUpdate<Derivative::DX, Derivative::NoOp, true, size_t>>,
  /* Hz */ FieldIntegratorNull<T>
>;

template<typename T>
// requires is_TM_fields<T, EmptyArray<typename T::value_t, T::dimension_t::value>>
using EMTM = TypeList<
  /* Ex */ FieldIntegratorNull<T>,
  /* Ey */ FieldIntegratorNull<T>,
  /* Ez */ FieldIntegrator2D<T, FieldUpdate<Derivative::DX, Derivative::DY, false, size_t, size_t>>,
  /* Hx */ FieldIntegrator2D<T, FieldUpdate<Derivative::NoOp, Derivative::DY, true, size_t, size_t>>,
  /* Hy */ FieldIntegrator2D<T, FieldUpdate<Derivative::DX, Derivative::NoOp, true>>,
  /* Hz */ FieldIntegratorNull<T>
>;

template<typename T>
// requires is_TM_fields<T, EmptyArray<typename T::value_t, T::dimension_t::value>>
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
  disable_field<EmptyArray1D<T>>, // Ex
  disable_field<EmptyArray1D<T>>, // Ey
  disable_field<EmptyArray1D<T>>, // Ez
  disable_field<EmptyArray1D<T>>, // Hx
  disable_field<EmptyArray1D<T>>, // Hy
  disable_field<EmptyArray1D<T>>  // Hz
>;

template<typename T>
using bcdata1d = BCData<
  disable_field<EmptyArray1D<T>>, // Ex
  disable_field<EmptyArray1D<T>>, // Ey
  enable_field<Array1D<T>>,       // Ez
  disable_field<EmptyArray1D<T>>, // Hx
  enable_field<Array1D<T>>,       // Hy
  disable_field<EmptyArray1D<T>>  // Hz
>;

template<typename T>
using bcdataTM = BCData<
  disable_field<EmptyArray2D<T>>, // Ex
  disable_field<EmptyArray2D<T>>, // Ey
  enable_field<Array2D<T>>,       // Ez
  enable_field<Array2D<T>>,       // Hx
  enable_field<Array2D<T>>,       // Hy
  disable_field<Array2D<T>>       // Hz
>;

template<typename T>
using bcdataTE = BCData<
  enable_field<Array2D<T>>,       // Ex
  enable_field<Array2D<T>>,       // Ey
  disable_field<EmptyArray2D<T>>, // Ez
  disable_field<EmptyArray2D<T>>, // Hx
  disable_field<EmptyArray2D<T>>, // Hy
  enable_field<Array2D<T>>        // Hz
>;

//=================== Boundary Conditions Definitions ========================
//=========================================================================
// template<typename EX, typename EY, typename EZ, typename HX, typename HY, typename HZ>
// // requires is_1D_fields<T, EmptyArray<typename T::value_t, T::dimension_t::value>>
// using BCNull = TypeList<
//   /* Ex */ NoneBC<EmptyArray1D<T>>,
//   /* Ey */ NoneBC<EmptyArray1D<T>>,
//   /* Ez */ NoneBC<EmptyArray1D<T>>,
//   /* Hx */ NoneBC<EmptyArray1D<T>>,
//   /* Hy */ NoneBC<EmptyArray1D<T>>,
//   /* Hz */ NoneBC<EmptyArray1D<T>>
// >;

#endif //EM_DEFINITIONS_H
