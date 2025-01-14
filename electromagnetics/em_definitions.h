//
// Created by cepheid on 10/17/24.
//

#ifndef EM_DEFINITIONS_H
#define EM_DEFINITIONS_H

#include "core/typelist.h"
#include "em_updates.h"
#include "em_data.h"
#include "em_solver.h"

namespace tf::electromagnetics
{
  //=================== EMData Definitions ========================
  //===============================================================
  template<typename T>
  using emdataNone = EMData<
    tf::types::EmptyArray3D<T>, // Ex
    tf::types::EmptyArray3D<T>, // Ey
    tf::types::EmptyArray3D<T>, // Ez
    tf::types::EmptyArray3D<T>, // Hx
    tf::types::EmptyArray3D<T>, // Hy
    tf::types::EmptyArray3D<T>  // Hz
  >;

  template<typename T>
  using emdata3D = EMData<
    tf::types::Array3D<T>, // Ex
    tf::types::Array3D<T>, // Ey
    tf::types::Array3D<T>, // Ez
    tf::types::Array3D<T>, // Hx
    tf::types::Array3D<T>, // Hy
    tf::types::Array3D<T>  // Hz
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
  using EM3D = TypeList<
    /* Ex */ FieldIntegrator3D<T, FieldUpdate<Derivative::DY, Derivative::DZ, true, size_t, size_t, size_t>>,
    /* Ey */ FieldIntegrator3D<T, FieldUpdate<Derivative::DZ, Derivative::DX, true, size_t, size_t, size_t>>,
    /* Ez */ FieldIntegrator3D<T, FieldUpdate<Derivative::DX, Derivative::DY, true, size_t, size_t, size_t>>,
    /* Hx */ FieldIntegrator3D<T, FieldUpdate<Derivative::DZ, Derivative::DY, false, size_t, size_t, size_t>>,
    /* Hy */ FieldIntegrator3D<T, FieldUpdate<Derivative::DX, Derivative::DZ, false, size_t, size_t, size_t>>,
    /* Hz */ FieldIntegrator3D<T, FieldUpdate<Derivative::DY, Derivative::DX, false, size_t, size_t, size_t>>
  >;
} // end namespace tf::electromagnetics
#endif //EM_DEFINITIONS_H
