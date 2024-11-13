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
  /* Ez */ PMLData<Array1D<T>, EMFace::X, S, false>,
  /* Hx */ NullData<EmptyArray1D<T>>,
  /* Hy */ PMLData<Array1D<T>, EMFace::X, S, true>,
  /* Hz */ NullData<EmptyArray1D<T>>
>;

template<typename T, EMSide S>
using PML_XFaceTM = FaceBCs<
  /* Ex */ NullData<EmptyArray2D<T>>,
  /* Ey */ NullData<EmptyArray2D<T>>,
  /* Ez */ PMLData<Array2D<T>, EMFace::X, S, false>,
  /* Hx */ NullData<EmptyArray2D<T>>,
  /* Hy */ PMLData<Array2D<T>, EMFace::X, S, true>,
  /* Hz */ NullData<EmptyArray2D<T>>
>;
template<typename T, EMSide S>
using PML_XFaceTE = FaceBCs<
  /* Ex */ NullData<EmptyArray2D<T>>,
  /* Ey */ PMLData<Array2D<T>, EMFace::X, S, false>,
  /* Ez */ NullData<EmptyArray2D<T>>,
  /* Hx */ NullData<EmptyArray2D<T>>,
  /* Hy */ NullData<EmptyArray2D<T>>,
  /* Hz */ PMLData<Array2D<T>, EMFace::X, S, true>
>;

template<typename T, EMSide S>
using PML_XFace3D = FaceBCs<
  /* Ex */ NullData<EmptyArray3D<T>>,
  /* Ey */ PMLData<Array3D<T>, EMFace::X, S, false>,
  /* Ez */ PMLData<Array3D<T>, EMFace::X, S, false>,
  /* Hx */ NullData<EmptyArray3D<T>>,
  /* Hy */ PMLData<Array3D<T>, EMFace::X, S, true>,
  /* Hz */ PMLData<Array3D<T>, EMFace::X, S, true>
>;


template<typename T, EMSide S>
using PML_YFaceTM = FaceBCs<
  /* Ex */ NullData<EmptyArray2D<T>>,
  /* Ey */ NullData<EmptyArray2D<T>>,
  /* Ez */ PMLData<Array2D<T>, EMFace::Y, S, false>,
  /* Hx */ PMLData<Array2D<T>, EMFace::Y, S, true>,
  /* Hy */ NullData<EmptyArray2D<T>>,
  /* Hz */ NullData<EmptyArray2D<T>>
>;

template<typename T, EMSide S>
using PML_YFaceTE = FaceBCs<
  /* Ex */ PMLData<Array2D<T>, EMFace::Y, S, false>,
  /* Ey */ NullData<EmptyArray2D<T>>,
  /* Ez */ NullData<EmptyArray2D<T>>,
  /* Hx */ NullData<EmptyArray2D<T>>,
  /* Hy */ NullData<EmptyArray2D<T>>,
  /* Hz */ PMLData<Array2D<T>, EMFace::Y, S, true>
>;

template<typename T, EMSide S>
using PML_YFace3D = FaceBCs<
  /* Ex */ PMLData<Array3D<T>, EMFace::Y, S, false>,
  /* Ey */ NullData<EmptyArray3D<T>>,
  /* Ez */ PMLData<Array3D<T>, EMFace::Y, S, false>,
  /* Hx */ PMLData<Array3D<T>, EMFace::Y, S, true>,
  /* Hy */ NullData<EmptyArray3D<T>>,
  /* Hz */ PMLData<Array3D<T>, EMFace::Y, S, true>
>;

template<typename T, EMSide S>
using PML_ZFace3D = FaceBCs<
  /* Ex */ PMLData<Array3D<T>, EMFace::Y, S, false>,
  /* Ey */ PMLData<Array3D<T>, EMFace::Y, S, false>,
  /* Ez */ NullData<EmptyArray3D<T>>,
  /* Hx */ PMLData<Array3D<T>, EMFace::Y, S, true>,
  /* Hy */ PMLData<Array3D<T>, EMFace::Y, S, true>,
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



#endif //BC_DEFINITIONS_H
