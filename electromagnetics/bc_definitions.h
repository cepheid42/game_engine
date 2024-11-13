//
// Created by cepheid on 11/11/24.
//

#ifndef BC_DEFINITIONS_H
#define BC_DEFINITIONS_H

#include <electromagnetics.h>

#include "core/typelist.h"
#include "bc_data.h"
#include "boundaries.h"

//=================== BCData Definitions ===================
//==========================================================
template<typename T, EMSide>
using Reflecting_Face = FaceBCs<
  /* Ex */ NullData<EmptyArray1D<T>>,
  /* Ey */ NullData<EmptyArray1D<T>>,
  /* Ez */ NullData<EmptyArray1D<T>>,
  /* Hx */ NullData<EmptyArray1D<T>>,
  /* Hy */ NullData<EmptyArray1D<T>>,
  /* Hz */ NullData<EmptyArray1D<T>>
>;

template<typename T, EMSide S>
using PML_XFace1D = FaceBCs<
  /* Ex */ NullData<EmptyArray1D<T>>,
  /* Ey */ NullData<EmptyArray1D<T>>,
  /* Ez */ PMLData<Array1D<T>, EMFace::X, S>,
  /* Hx */ NullData<EmptyArray1D<T>>,
  /* Hy */ PMLData<Array1D<T>, EMFace::X, S>,
  /* Hz */ NullData<EmptyArray1D<T>>
>;

template<typename T, EMSide S>
using PML_XFaceTM = FaceBCs<
  /* Ex */ NullData<EmptyArray2D<T>>,
  /* Ey */ NullData<EmptyArray2D<T>>,
  /* Ez */ PMLData<Array2D<T>, EMFace::X, S>,
  /* Hx */ NullData<EmptyArray2D<T>>,
  /* Hy */ PMLData<Array2D<T>, EMFace::X, S>,
  /* Hz */ NullData<EmptyArray2D<T>>
>;
template<typename T, EMSide S>
using PML_XFaceTE = FaceBCs<
  /* Ex */ NullData<EmptyArray2D<T>>,
  /* Ey */ PMLData<Array2D<T>, EMFace::X, S>,
  /* Ez */ NullData<EmptyArray2D<T>>,
  /* Hx */ NullData<EmptyArray2D<T>>,
  /* Hy */ NullData<EmptyArray2D<T>>,
  /* Hz */ PMLData<Array2D<T>, EMFace::X, S>
>;

template<typename T, EMSide S>
using PML_XFace3D = FaceBCs<
  /* Ex */ NullData<EmptyArray3D<T>>,
  /* Ey */ PMLData<Array3D<T>, EMFace::X, S>,
  /* Ez */ PMLData<Array3D<T>, EMFace::X, S>,
  /* Hx */ NullData<EmptyArray3D<T>>,
  /* Hy */ PMLData<Array3D<T>, EMFace::X, S>,
  /* Hz */ PMLData<Array3D<T>, EMFace::X, S>
>;


template<typename T, EMSide S>
using PML_YFaceTM = FaceBCs<
  /* Ex */ NullData<EmptyArray2D<T>>,
  /* Ey */ NullData<EmptyArray2D<T>>,
  /* Ez */ PMLData<Array2D<T>, EMFace::Y, S>,
  /* Hx */ PMLData<Array2D<T>, EMFace::Y, S>,
  /* Hy */ NullData<EmptyArray2D<T>>,
  /* Hz */ NullData<EmptyArray2D<T>>
>;

template<typename T, EMSide S>
using PML_YFaceTE = FaceBCs<
  /* Ex */ PMLData<Array2D<T>, EMFace::Y, S>,
  /* Ey */ NullData<EmptyArray2D<T>>,
  /* Ez */ NullData<EmptyArray2D<T>>,
  /* Hx */ NullData<EmptyArray2D<T>>,
  /* Hy */ NullData<EmptyArray2D<T>>,
  /* Hz */ PMLData<Array2D<T>, EMFace::Y, S>
>;

template<typename T, EMSide S>
using PML_YFace3D = FaceBCs<
  /* Ex */ PMLData<Array3D<T>, EMFace::Y, S>,
  /* Ey */ NullData<EmptyArray3D<T>>,
  /* Ez */ PMLData<Array3D<T>, EMFace::Y, S>,
  /* Hx */ PMLData<Array3D<T>, EMFace::Y, S>,
  /* Hy */ NullData<EmptyArray3D<T>>,
  /* Hz */ PMLData<Array3D<T>, EMFace::Y, S>
>;

template<typename T, EMSide S>
using PML_ZFace3D = FaceBCs<
  /* Ex */ PMLData<Array3D<T>, EMFace::Y, S>,
  /* Ey */ PMLData<Array3D<T>, EMFace::Y, S>,
  /* Ez */ NullData<EmptyArray3D<T>>,
  /* Hx */ PMLData<Array3D<T>, EMFace::Y, S>,
  /* Hy */ PMLData<Array3D<T>, EMFace::Y, S>,
  /* Hz */ NullData<EmptyArray3D<T>>
>;

template<typename T, EMFace F, EMSide S>
using Periodic_Face1D = FaceBCs<
  /* Ex */ NullData<EmptyArray1D<T>>,
  /* Ey */ NullData<EmptyArray1D<T>>,
  /* Ez */ NullData<EmptyArray1D<T>>,
  /* Hx */ NullData<EmptyArray1D<T>>,
  /* Hy */ PeriodicData<Array1D<T>, F, S>,
  /* Hz */ NullData<EmptyArray1D<T>>
>;

template<typename T, EMFace F, EMSide S>
using Periodic_FaceTM = FaceBCs<
  /* Ex */ NullData<EmptyArray2D<T>>,
  /* Ey */ NullData<EmptyArray2D<T>>,
  /* Ez */ NullData<EmptyArray2D<T>>,
  /* Hx */ PeriodicData<Array2D<T>, F, S>,
  /* Hy */ PeriodicData<Array2D<T>, F, S>,
  /* Hz */ NullData<EmptyArray2D<T>>
>;

template<typename T, EMFace F, EMSide S>
using Periodic_FaceTE = FaceBCs<
  /* Ex */ NullData<EmptyArray2D<T>>,
  /* Ey */ NullData<EmptyArray2D<T>>,
  /* Ez */ NullData<EmptyArray2D<T>>,
  /* Hx */ NullData<EmptyArray2D<T>>,
  /* Hy */ NullData<EmptyArray2D<T>>,
  /* Hz */ PeriodicData<Array2D<T>, F, S>
>;

template<typename T, EMFace F, EMSide S>
using Periodic_Face3D = FaceBCs<
  /* Ex */ NullData<EmptyArray3D<T>>,
  /* Ey */ NullData<EmptyArray3D<T>>,
  /* Ez */ NullData<EmptyArray3D<T>>,
  /* Hx */ PeriodicData<Array3D<T>, F, S>,
  /* Hy */ PeriodicData<Array3D<T>, F, S>,
  /* Hz */ PeriodicData<Array3D<T>, F, S>
>;

//=================== Boundary Condition Definitions ===================
//======================================================================
template<typename T>
using BCReflecting = EMBoundary<
  EMFace::X,
  /* Ex */ BCIntegratorNull<T>,
  /* Ey */ BCIntegratorNull<T>,
  /* Ez */ BCIntegratorNull<T>,
  /* Hx */ BCIntegratorNull<T>,
  /* Hy */ BCIntegratorNull<T>,
  /* Hz */ BCIntegratorNull<T>
>;

template<typename T>
using BCPeriodic1D = EMBoundary<
  EMFace::X,
  /* Ex */ BCIntegratorNull<T>,
  /* Ey */ BCIntegratorNull<T>,
  /* Ez */ BCIntegratorNull<T>,
  /* Hx */ BCIntegratorNull<T>,
  /* Hy */ BCIntegrator<T, PeriodicBC<size_t>>,
  /* Hz */ BCIntegratorNull<T>
>;

template<typename T>
using BCPeriodicTM = EMBoundary<
  EMFace::X,
  /* Ex */ BCIntegratorNull<T>,
  /* Ey */ BCIntegratorNull<T>,
  /* Ez */ BCIntegrator<T, PeriodicBC<size_t, size_t>>,
  /* Hx */ BCIntegrator<T, PeriodicBC<size_t, size_t>>,
  /* Hy */ BCIntegrator<T, PeriodicBC<size_t, size_t>>,
  /* Hz */ BCIntegratorNull<T>
>;

template<typename T>
using BCPeriodicTE = EMBoundary<
  EMFace::X,
  /* Ex */ BCIntegrator<T, PeriodicBC<size_t, size_t>>,
  /* Ey */ BCIntegrator<T, PeriodicBC<size_t, size_t>>,
  /* Ez */ BCIntegratorNull<T>,
  /* Hx */ BCIntegratorNull<T>,
  /* Hy */ BCIntegratorNull<T>,
  /* Hz */ BCIntegrator<T, PeriodicBC<size_t, size_t>>
>;

template<typename T>
using BCPeriodic3D = EMBoundary<
  EMFace::X,
  /* Ex */ BCIntegrator<T, PeriodicBC<size_t, size_t, size_t>>,
  /* Ey */ BCIntegrator<T, PeriodicBC<size_t, size_t, size_t>>,
  /* Ez */ BCIntegrator<T, PeriodicBC<size_t, size_t, size_t>>,
  /* Hx */ BCIntegrator<T, PeriodicBC<size_t, size_t, size_t>>,
  /* Hy */ BCIntegrator<T, PeriodicBC<size_t, size_t, size_t>>,
  /* Hz */ BCIntegrator<T, PeriodicBC<size_t, size_t, size_t>>
>;

template<typename T>
using BCPml1D = EMBoundary<
  EMFace::X,
  /* Ex */ BCIntegratorNull<T>,
  /* Ey */ BCIntegratorNull<T>,
  /* Ez */ BCIntegrator<T, PmlBC<Derivative::DX, Derivative::NoOp, false, size_t>>,
  /* Hx */ BCIntegratorNull<T>,
  /* Hy */ BCIntegrator<T, PmlBC<Derivative::NoOp, Derivative::DX, true, size_t>>,
  /* Hz */ BCIntegratorNull<T>
>;

template<typename T>
using BCPmlTM_XFace = EMBoundary<
  EMFace::X,
  /* Ex */ BCIntegratorNull<T>,
  /* Ey */ BCIntegratorNull<T>,
  /* Ez */ BCIntegrator<T, PmlBC<Derivative::DX, Derivative::NoOp, false, size_t, size_t>>,
  /* Hx */ BCIntegratorNull<T>,
  /* Hy */ BCIntegrator<T, PmlBC<Derivative::NoOp, Derivative::DX, true, size_t, size_t>>,
  /* Hz */ BCIntegratorNull<T>
>;

template<typename T>
using BCPmlTM_YFace = EMBoundary<
  EMFace::Y,
  /* Ex */ BCIntegratorNull<T>,
  /* Ey */ BCIntegratorNull<T>,
  /* Ez */ BCIntegrator<T, PmlBC<Derivative::NoOp, Derivative::DY, false, size_t, size_t>>,
  /* Hx */ BCIntegrator<T, PmlBC<Derivative::DY, Derivative::NoOp, true, size_t, size_t>>,
  /* Hy */ BCIntegratorNull<T>,
  /* Hz */ BCIntegratorNull<T>
>;

template<typename T>
using BCPmlTE_XFace = EMBoundary<
  EMFace::X,
  /* Ex */ BCIntegratorNull<T>,
  /* Ey */ BCIntegrator<T, PmlBC<Derivative::NoOp, Derivative::DX, false, size_t, size_t>>,
  /* Ez */ BCIntegratorNull<T>,
  /* Hx */ BCIntegratorNull<T>,
  /* Hy */ BCIntegratorNull<T>,
  /* Hz */ BCIntegrator<T, PmlBC<Derivative::DX, Derivative::NoOp, true, size_t, size_t>>
>;

template<typename T>
using BCPmlTE_YFace = EMBoundary<
  EMFace::Y,
  /* Ex */ BCIntegrator<T, PmlBC<Derivative::DY, Derivative::NoOp, false, size_t, size_t>>,
  /* Ey */ BCIntegratorNull<T>,
  /* Ez */ BCIntegratorNull<T>,
  /* Hx */ BCIntegratorNull<T>,
  /* Hy */ BCIntegratorNull<T>,
  /* Hz */ BCIntegrator<T, PmlBC<Derivative::NoOp, Derivative::DY, true, size_t, size_t>>
>;

template<typename T>
using BCPml3D_XFace = EMBoundary<
  EMFace::X,
  /* Ex */ BCIntegratorNull<T>,
  /* Ey */ BCIntegrator<T, PmlBC<Derivative::NoOp, Derivative::DX, false, size_t, size_t, size_t>>,
  /* Ez */ BCIntegrator<T, PmlBC<Derivative::DX, Derivative::NoOp, false, size_t, size_t, size_t>>,
  /* Hx */ BCIntegratorNull<T>,
  /* Hy */ BCIntegrator<T, PmlBC<Derivative::NoOp, Derivative::DX, true, size_t, size_t, size_t>>,
  /* Hz */ BCIntegrator<T, PmlBC<Derivative::DX, Derivative::NoOp, true, size_t, size_t, size_t>>
>;

template<typename T>
using BCPml3D_YFace = EMBoundary<
  EMFace::Y,
  /* Ex */ BCIntegrator<T, PmlBC<Derivative::DY, Derivative::NoOp, false, size_t, size_t, size_t>>,
  /* Ey */ BCIntegratorNull<T>,
  /* Ez */ BCIntegrator<T, PmlBC<Derivative::NoOp, Derivative::DY, false, size_t, size_t, size_t>>,
  /* Hx */ BCIntegrator<T, PmlBC<Derivative::DY, Derivative::NoOp, true, size_t, size_t, size_t>>,
  /* Hy */ BCIntegratorNull<T>,
  /* Hz */ BCIntegrator<T, PmlBC<Derivative::NoOp, Derivative::DY, true, size_t, size_t, size_t>>
>;

template<typename T>
using BCPml3D_ZFace = EMBoundary<
  EMFace::Z,
  /* Ex */ BCIntegrator<T, PmlBC<Derivative::NoOp, Derivative::DZ, false, size_t, size_t, size_t>>,
  /* Ey */ BCIntegrator<T, PmlBC<Derivative::DZ, Derivative::NoOp, false, size_t, size_t, size_t>>,
  /* Ez */ BCIntegratorNull<T>,
  /* Hx */ BCIntegrator<T, PmlBC<Derivative::NoOp, Derivative::DZ, true, size_t, size_t, size_t>>,
  /* Hy */ BCIntegrator<T, PmlBC<Derivative::DZ, Derivative::NoOp, true, size_t, size_t, size_t>>,
  /* Hz */ BCIntegratorNull<T>
>;


#endif //BC_DEFINITIONS_H
