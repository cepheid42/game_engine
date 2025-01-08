//
// Created by cepheid on 10/17/24.
//

#ifndef EM_DEFINITIONS_H
#define EM_DEFINITIONS_H

#include "core/typelist.h"
#include "em_updates.h"
#include "em_data.h"
#include "em_solver.h"

// todo: Does this need to have namespaces? Are these usings polluting upstream?
using tf::electromagnetics::traits::Derivative;
using tf::electromagnetics::EMData;
using tf::electromagnetics::FieldUpdate;
using tf::electromagnetics::FieldIntegratorNull;
using tf::electromagnetics::FieldIntegrator3D;
using tf::types::Array3D;
using tf::types::EmptyArray3D;

//=================== EMData Definitions ========================
//===============================================================
template<typename T>
using emdataNone = EMData<
  EmptyArray3D<T>, // Ex
  EmptyArray3D<T>, // Ey
  EmptyArray3D<T>, // Ez
  EmptyArray3D<T>, // Hx
  EmptyArray3D<T>, // Hy
  EmptyArray3D<T>  // Hz
>;

template<typename T>
using emdata1D = EMData<
  EmptyArray3D<T>, // Ex
  EmptyArray3D<T>, // Ey
  Array3D<T>,      // Ez
  EmptyArray3D<T>, // Hx
  Array3D<T>,      // Hy
  EmptyArray3D<T>  // Hz
>;

template<typename T>
using emdataTMz = EMData<
  EmptyArray3D<T>, // Ex
  EmptyArray3D<T>, // Ey
  Array3D<T>,      // Ez
  Array3D<T>,      // Hx
  Array3D<T>,      // Hy
  EmptyArray3D<T>  // Hz
>;

template<typename T>
using emdataTEz = EMData<
  Array3D<T>,      // Ex
  Array3D<T>,      // Ey
  EmptyArray3D<T>, // Ez
  EmptyArray3D<T>, // Hx
  EmptyArray3D<T>, // Hy
  Array3D<T>       // Hz
>;

template<typename T>
using emdata3D = EMData<
  Array3D<T>, // Ex
  Array3D<T>, // Ey
  Array3D<T>, // Ez
  Array3D<T>, // Hx
  Array3D<T>, // Hy
  Array3D<T>  // Hz
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
  /* Ez */ FieldIntegrator3D<T, FieldUpdate<Derivative::DX, Derivative::NoOp, false, size_t, size_t, size_t>>,
  /* Hx */ FieldIntegratorNull<T>,
  /* Hy */ FieldIntegrator3D<T, FieldUpdate<Derivative::DX, Derivative::NoOp, true, size_t, size_t, size_t>>,
  /* Hz */ FieldIntegratorNull<T>
>;

template<typename T>
using EMTMz = TypeList<
  /* Ex */ FieldIntegratorNull<T>,
  /* Ey */ FieldIntegratorNull<T>,
  /* Ez */ FieldIntegrator3D<T, FieldUpdate<Derivative::DX, Derivative::DY, false, size_t, size_t, size_t>>,
  /* Hx */ FieldIntegrator3D<T, FieldUpdate<Derivative::NoOp, Derivative::DY, true, size_t, size_t, size_t>>,
  /* Hy */ FieldIntegrator3D<T, FieldUpdate<Derivative::DX, Derivative::NoOp, true, size_t, size_t, size_t>>,
  /* Hz */ FieldIntegratorNull<T>
>;

template<typename T>
using EMTEz = TypeList<
  /* Ex */ FieldIntegrator3D<T, FieldUpdate<Derivative::DY, Derivative::NoOp, false, size_t, size_t, size_t>>,
  /* Ey */ FieldIntegrator3D<T, FieldUpdate<Derivative::NoOp, Derivative::DX, false, size_t, size_t, size_t>>,
  /* Ez */ FieldIntegratorNull<T>,
  /* Hx */ FieldIntegratorNull<T>,
  /* Hy */ FieldIntegratorNull<T>,
  /* Hz */ FieldIntegrator3D<T, FieldUpdate<Derivative::DY, Derivative::DX, true, size_t, size_t, size_t>>
>;

template<typename T>
using EM3D = TypeList<
  /* Ex */ FieldIntegrator3D<T, FieldUpdate<Derivative::DY, Derivative::DZ, false, size_t, size_t, size_t>>,
  /* Ey */ FieldIntegrator3D<T, FieldUpdate<Derivative::DZ, Derivative::DX, false, size_t, size_t, size_t>>,
  /* Ez */ FieldIntegrator3D<T, FieldUpdate<Derivative::DX, Derivative::DY, false, size_t, size_t, size_t>>,
  /* Hx */ FieldIntegrator3D<T, FieldUpdate<Derivative::DZ, Derivative::DY, true, size_t, size_t, size_t>>,
  /* Hy */ FieldIntegrator3D<T, FieldUpdate<Derivative::DX, Derivative::DZ, true, size_t, size_t, size_t>>,
  /* Hz */ FieldIntegrator3D<T, FieldUpdate<Derivative::DY, Derivative::DX, true, size_t, size_t, size_t>>
>;

#endif //EM_DEFINITIONS_H
