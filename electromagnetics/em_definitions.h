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
using tf::electromagnetics::FieldIntegrator1D;
using tf::electromagnetics::FieldIntegrator2D;
using tf::electromagnetics::FieldIntegrator3D;
using tf::types::Array1D;
using tf::types::Array2D;
using tf::types::Array3D;
using tf::types::EmptyArray1D;
using tf::types::EmptyArray2D;
using tf::types::EmptyArray3D;

//=================== EMData Definitions ========================
//===============================================================
template<typename T>
using emdataNone = EMData<
  EmptyArray1D<T>, // Ex
  EmptyArray1D<T>, // Ey
  EmptyArray1D<T>, // Ez
  EmptyArray1D<T>, // Hx
  EmptyArray1D<T>, // Hy
  EmptyArray1D<T>  // Hz
>;

template<typename T>
using emdata1D = EMData<
  EmptyArray1D<T>, // Ex
  EmptyArray1D<T>, // Ey
  Array1D<T>,      // Ez
  EmptyArray1D<T>, // Hx
  Array1D<T>,      // Hy
  EmptyArray1D<T>  // Hz
>;

template<typename T>
using emdataTM = EMData<
  EmptyArray2D<T>, // Ex
  EmptyArray2D<T>, // Ey
  Array2D<T>,      // Ez
  Array2D<T>,      // Hx
  Array2D<T>,      // Hy
  EmptyArray2D<T>  // Hz
>;

template<typename T>
using emdataTE = EMData<
  Array2D<T>,      // Ex
  Array2D<T>,      // Ey
  EmptyArray2D<T>, // Ez
  EmptyArray2D<T>, // Hx
  EmptyArray2D<T>, // Hy
  Array2D<T>       // Hz
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
