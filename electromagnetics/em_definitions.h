//
// Created by cepheid on 10/17/24.
//

#ifndef EM_DEFINITIONS_H
#define EM_DEFINITIONS_H

#include "core/typelist.h"
#include "em_data.h"
#include "em_solver.h"
#include "bc_data.h"
#include "boundaries.h"
#include "offsets.h"

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
  Array1D<T>,       // Ez
  EmptyArray1D<T>, // Hx
  Array1D<T>,       // Hy
  EmptyArray1D<T>  // Hz
>;

template<typename T>
using emdataTM = EMData<
  EmptyArray2D<T>, // Ex
  EmptyArray2D<T>, // Ey
  Array2D<T>,       // Ez
  Array2D<T>,       // Hx
  Array2D<T>,       // Hy
  EmptyArray2D<T>  // Hz
>;

template<typename T>
using emdataTE = EMData<
  Array2D<T>,       // Ex
  Array2D<T>,       // Ey
  EmptyArray2D<T>, // Ez
  EmptyArray2D<T>, // Hx
  EmptyArray2D<T>, // Hy
  Array2D<T>        // Hz
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

//=================== Boundary Condition Definitions ========================
//===========================================================================

template<typename Array, bool Forward, typename... IDXS>
using EyHy_XFace_PML = PmlBC<
  Array,
  curl<Derivative::DZ, Forward, IDXS...>,
  curl<Derivative::NoOp, Forward, IDXS...>,
  dPML,
  IDXS...
>;

template<typename Array, bool Forward, typename... IDXS>
using EzHz_XFace_PML = PmlBC<
  Array,
  curl<Derivative::NoOp, Forward, IDXS...>,
  curl<Derivative::DY, Forward, IDXS...>,
  dPML,
  IDXS...
>;

template<typename Array, bool Forward, typename... IDXS>
using ExHx_YFace_PML = PmlBC<
  Array,
  curl<Derivative::NoOp, Forward, IDXS...>,
  curl<Derivative::DZ, Forward, IDXS...>,
  dPML,
  IDXS...
>;

template<typename Array, bool Forward, typename... IDXS>
using EzHz_YFace_PML = PmlBC<
  Array,
  curl<Derivative::DX, Forward, IDXS...>,
  curl<Derivative::NoOp, Forward, IDXS...>,
  dPML,
  IDXS...
>;

template<typename Array, bool Forward, typename... IDXS>
using ExHx_ZFace_PML = PmlBC<
  Array,
  curl<Derivative::DY, Forward, IDXS...>,
  curl<Derivative::NoOp, Forward, IDXS...>,
  dPML,
  IDXS...
>;

template<typename Array, bool Forward, typename... IDXS>
using EyHy_ZFace_PML = PmlBC<
  Array,
  curl<Derivative::NoOp, Forward, IDXS...>,
  curl<Derivative::DX, Forward, IDXS...>,
  dPML,
  IDXS...
>;

template<typename Array, bool Forward, typename... IDXS>
using XFaceTL = TypeList<
  EyHy_XFace_PML<Array, Forward, IDXS...>,
  EzHz_XFace_PML<Array, Forward, IDXS...>,
  PeriodicBC<Array, EMFace::X, nHalo, IDXS...>,
  ReflectingBC<Array>
>;

template<typename Array, bool Forward, typename... IDXS>
using YFaceTL = TypeList<
  ExHx_YFace_PML<Array, Forward, IDXS...>,
  EzHz_YFace_PML<Array, Forward, IDXS...>,
  PeriodicBC<Array, EMFace::Y, nHalo, IDXS...>,
  ReflectingBC<Array>
>;

template<typename Array, bool Forward, typename... IDXS>
using ZFaceTL = TypeList<
  ExHx_ZFace_PML<Array, Forward, IDXS...>,
  EyHy_ZFace_PML<Array, Forward, IDXS...>,
  PeriodicBC<Array, EMFace::Z, nHalo, IDXS...>,
  ReflectingBC<Array>
>;


template<size_t I, typename T>
using X0Face_1D = BCIntegrator1D<TypeListAt<I, XFaceTL<Array1D<T>, false, size_t>>>;

template<size_t I, typename T>
using X1Face_1D = BCIntegrator1D<TypeListAt<I, XFaceTL<Array1D<T>, false, size_t>>>;

template<typename T>
using BoundariesTL = TypeList<
  X0Face_1D<1, T>, // Ez X0
  X1Face_1D<1, T>  // Ez X1
>;

#endif //EM_DEFINITIONS_H
