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
  /* Ez */ FieldIntegrator1D<T, Derivative::DX, Derivative::NoOp, false>,
  /* Hx */ FieldIntegratorNull<T>,
  /* Hy */ FieldIntegrator1D<T, Derivative::DX, Derivative::NoOp, true>,
  /* Hz */ FieldIntegratorNull<T>
>;

template<typename T>
// requires is_TM_fields<T, EmptyArray<typename T::value_t, T::dimension_t::value>>
using EMTM = TypeList<
  /* Ex */ FieldIntegratorNull<T>,
  /* Ey */ FieldIntegratorNull<T>,
  /* Ez */ FieldIntegrator2D<T, Derivative::DX, Derivative::DY, false>,
  /* Hx */ FieldIntegrator2D<T, Derivative::NoOp, Derivative::DY, true>,
  /* Hy */ FieldIntegrator2D<T, Derivative::DX, Derivative::NoOp, true>,
  /* Hz */ FieldIntegratorNull<T>
>;

#endif //EM_DEFINITIONS_H
