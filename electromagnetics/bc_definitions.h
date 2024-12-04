//
// Created by cepheid on 11/11/24.
//

#ifndef BC_DEFINITIONS_H
#define BC_DEFINITIONS_H

#include "bc_data.h"
#include "boundaries.h"

//=================== BCData Definitions ===================
//==========================================================
// -------------------------------------------
// Reflecting Boundary
template<typename T, EMSide>
using ReflectingFace = FaceBCs<
  /* Ex */ NullData<EmptyArray1D<T>>,
  /* Ey */ NullData<EmptyArray1D<T>>,
  /* Ez */ NullData<EmptyArray1D<T>>,
  /* Hx */ NullData<EmptyArray1D<T>>,
  /* Hy */ NullData<EmptyArray1D<T>>,
  /* Hz */ NullData<EmptyArray1D<T>>
>;

// -------------------------------------------
// 1D Boundary
template<typename T, EMFace F, EMSide S>
using PeriodicFace1D = FaceBCs<
  /* Ex */ NullData<EmptyArray1D<T>>,
  /* Ey */ NullData<EmptyArray1D<T>>,
  /* Ez */ NullData<EmptyArray1D<T>>,
  /* Hx */ NullData<EmptyArray1D<T>>,
  /* Hy */ PeriodicData<Array1D<T>, F, S>,
  /* Hz */ NullData<EmptyArray1D<T>>
>;

// -------------------------------------------
// TM Boundary
template<typename T, EMFace F, EMSide S>
using PeriodicFaceTM = FaceBCs<
  /* Ex */ NullData<EmptyArray2D<T>>,
  /* Ey */ NullData<EmptyArray2D<T>>,
  /* Ez */ NullData<EmptyArray2D<T>>,
  /* Hx */ PeriodicData<Array2D<T>, F, S>,
  /* Hy */ PeriodicData<Array2D<T>, F, S>,
  /* Hz */ NullData<EmptyArray2D<T>>
>;

// -------------------------------------------
// TE Boundary
template<typename T, EMFace F, EMSide S>
using PeriodicFaceTE = FaceBCs<
  /* Ex */ NullData<EmptyArray2D<T>>,
  /* Ey */ NullData<EmptyArray2D<T>>,
  /* Ez */ NullData<EmptyArray2D<T>>,
  /* Hx */ NullData<EmptyArray2D<T>>,
  /* Hy */ NullData<EmptyArray2D<T>>,
  /* Hz */ PeriodicData<Array2D<T>, F, S>
>;

// -------------------------------------------
// 3D Periodic Boundary
template<typename T, EMFace F, EMSide S>
using PeriodicFace3D = FaceBCs<
  /* Ex */ NullData<EmptyArray3D<T>>,
  /* Ey */ NullData<EmptyArray3D<T>>,
  /* Ez */ NullData<EmptyArray3D<T>>,
  /* Hx */ PeriodicData<Array3D<T>, F, S>,
  /* Hy */ PeriodicData<Array3D<T>, F, S>,
  /* Hz */ PeriodicData<Array3D<T>, F, S>
>;

// -------------------------------------------
// 1D PML Boundary
template<typename T, EMFace F, EMSide S>
using PmlFace1D = FaceBCs<
  /* Ex */ NullData<EmptyArray1D<T>>,
  /* Ey */ NullData<EmptyArray1D<T>>,
  /* Ez */ PMLData<Array1D<T>, F, S, false>,
  /* Hx */ NullData<EmptyArray1D<T>>,
  /* Hy */ PMLData<Array1D<T>, F, S, true>,
  /* Hz */ NullData<EmptyArray1D<T>>
>;

// -------------------------------------------
// TM PML
template<typename T, EMFace F, EMSide S>
struct PmlFaceTMImpl {
  using type = FaceBCs<
  /* Ex */ NullData<EmptyArray2D<T>>,
  /* Ey */ NullData<EmptyArray2D<T>>,
  /* Ez */ NullData<EmptyArray2D<T>>,
  /* Hx */ NullData<EmptyArray2D<T>>,
  /* Hy */ NullData<EmptyArray2D<T>>,
  /* Hz */ NullData<EmptyArray2D<T>>
>;
};

// TM X-Face
template<typename T, EMFace F, EMSide S>
requires (F == EMFace::X)
struct PmlFaceTMImpl<T, F, S> {
  using type = FaceBCs<
    /* Ex */ NullData<EmptyArray2D<T>>,
    /* Ey */ NullData<EmptyArray2D<T>>,
    /* Ez */ PMLData<Array2D<T>, F, S, false>,
    /* Hx */ NullData<EmptyArray2D<T>>,
    /* Hy */ PMLData<Array2D<T>, F, S, true>,
    /* Hz */ NullData<EmptyArray2D<T>>
  >;
};

// TM Y-Face
template<typename T, EMFace F, EMSide S>
requires (F == EMFace::Y)
struct PmlFaceTMImpl<T, F, S> {
  using type = FaceBCs<
    /* Ex */ NullData<EmptyArray2D<T>>,
    /* Ey */ NullData<EmptyArray2D<T>>,
    /* Ez */ PMLData<Array2D<T>, F, S, false>,
    /* Hx */ PMLData<Array2D<T>, F, S, true>,
    /* Hy */ NullData<EmptyArray2D<T>>,
    /* Hz */ NullData<EmptyArray2D<T>>
  >;
};

// Top-level alias
template<typename T, EMFace F, EMSide S>
using PmlFaceTM = typename PmlFaceTMImpl<T, F, S>::type;


// -------------------------------------------
// TE PML
template<typename T, EMFace F, EMSide S>
struct PmlFaceTEImpl;

// TE X-Face
template<typename T, EMFace F, EMSide S>
requires (F == EMFace::X)
struct PmlFaceTEImpl<T, F, S> {
  using type = FaceBCs<
    /* Ex */ NullData<EmptyArray2D<T>>,
    /* Ey */ PMLData<Array2D<T>, F, S, false>,
    /* Ez */ NullData<EmptyArray2D<T>>,
    /* Hx */ NullData<EmptyArray2D<T>>,
    /* Hy */ NullData<EmptyArray2D<T>>,
    /* Hz */ PMLData<Array2D<T>, F, S, true>
  >;
};

// TE Y-Face
template<typename T, EMFace F, EMSide S>
requires (F == EMFace::Y)
struct PmlFaceTEImpl<T, F, S> {
  using type = FaceBCs<
    /* Ex */ PMLData<Array2D<T>, F, S, false>,
    /* Ey */ NullData<EmptyArray2D<T>>,
    /* Ez */ NullData<EmptyArray2D<T>>,
    /* Hx */ NullData<EmptyArray2D<T>>,
    /* Hy */ NullData<EmptyArray2D<T>>,
    /* Hz */ PMLData<Array2D<T>, F, S, true>
  >;
};

// Top-level alias
template<typename T, EMFace F, EMSide S>
using PmlFaceTE = typename PmlFaceTMImpl<T, F, S>::type;

// -------------------------------------------
// 3D Pml
template<typename T, EMFace F, EMSide S>
struct PmlFace3DImpl;

// 3D X-Face
template<typename T, EMFace F, EMSide S>
requires (F == EMFace::X)
struct PmlFace3DImpl<T, F, S> {
  using type = FaceBCs<
    /* Ex */ NullData<EmptyArray3D<T>>,
    /* Ey */ PMLData<Array3D<T>, F, S, false>,
    /* Ez */ PMLData<Array3D<T>, F, S, false>,
    /* Hx */ NullData<EmptyArray3D<T>>,
    /* Hy */ PMLData<Array3D<T>, F, S, true>,
    /* Hz */ PMLData<Array3D<T>, F, S, true>
  >;
};

// 3D Y-Face
template<typename T, EMFace F, EMSide S>
requires (F == EMFace::Y)
struct PmlFace3DImpl<T, F, S> {
  using type = FaceBCs<
    /* Ex */ PMLData<Array3D<T>, F, S, false>,
    /* Ey */ NullData<EmptyArray3D<T>>,
    /* Ez */ PMLData<Array3D<T>, F, S, false>,
    /* Hx */ PMLData<Array3D<T>, F, S, true>,
    /* Hy */ NullData<EmptyArray3D<T>>,
    /* Hz */ PMLData<Array3D<T>, F, S, true>
  >;
};

// 3D Z-Face
template<typename T, EMFace F, EMSide S>
requires (F == EMFace::Z)
struct PmlFace3DImpl<T, F, S> {
  using type = FaceBCs<
    /* Ex */ PMLData<Array3D<T>, F, S, false>,
    /* Ey */ PMLData<Array3D<T>, F, S, false>,
    /* Ez */ NullData<EmptyArray3D<T>>,
    /* Hx */ PMLData<Array3D<T>, F, S, true>,
    /* Hy */ PMLData<Array3D<T>, F, S, true>,
    /* Hz */ NullData<EmptyArray3D<T>>
  >;
};

// Top-level alias
template<typename T, EMFace F, EMSide S>
using PmlFace3D = typename PmlFace3DImpl<T, F, S>::type;

//=================== Boundary Condition Definitions ===================
//======================================================================
using BCNull = TypeList<
  /* Ex */ ReflectingBCUpdate,
  /* Ey */ ReflectingBCUpdate,
  /* Ez */ ReflectingBCUpdate,
  /* Hx */ ReflectingBCUpdate,
  /* Hy */ ReflectingBCUpdate,
  /* Hz */ ReflectingBCUpdate
>;

using Periodic1D = TypeList<
  /* Ex */ Periodic1DUpdate,
  /* Ey */ Periodic1DUpdate,
  /* Ez */ Periodic1DUpdate,
  /* Hx */ Periodic1DUpdate,
  /* Hy */ Periodic1DUpdate,
  /* Hz */ Periodic1DUpdate
>;

template<EMFace F, EMSide S>
using Periodic2D = TypeList<
  /* Ex */ Periodic2DUpdate<F, S>,
  /* Ey */ Periodic2DUpdate<F, S>,
  /* Ez */ Periodic2DUpdate<F, S>,
  /* Hx */ Periodic2DUpdate<F, S>,
  /* Hy */ Periodic2DUpdate<F, S>,
  /* Hz */ Periodic2DUpdate<F, S>
>;


template<EMFace F, EMSide S>
using Periodic3D = TypeList<
  /* Ex */ Periodic3DUpdate<F, S>,
  /* Ey */ Periodic3DUpdate<F, S>,
  /* Ez */ Periodic3DUpdate<F, S>,
  /* Hx */ Periodic3DUpdate<F, S>,
  /* Hy */ Periodic3DUpdate<F, S>,
  /* Hz */ Periodic3DUpdate<F, S>
>;

template<EMSide S>
using PmlX1D = TypeList<
  /* Ex */ ReflectingBCUpdate,
  /* Ey */ ReflectingBCUpdate,
  /* Ez */ Pml1DUpdate<S>,
  /* Hx */ ReflectingBCUpdate,
  /* Hy */ Pml1DUpdate<S>,
  /* Hz */ ReflectingBCUpdate
>;

template<EMSide S>
using PmlXTM = TypeList<
  /* Ex */ ReflectingBCUpdate,
  /* Ey */ ReflectingBCUpdate,
  /* Ez */ Pml3DUpdate<EMFace::X, S, false>,
  /* Hx */ ReflectingBCUpdate,
  /* Hy */ Pml3DUpdate<EMFace::X, S, false>,
  /* Hz */ ReflectingBCUpdate
>;

template<EMSide S>
using PmlYTM = TypeList<
  /* Ex */ ReflectingBCUpdate,
  /* Ey */ ReflectingBCUpdate,
  /* Ez */ Pml3DUpdate<EMFace::Y, S, true>,
  /* Hx */ Pml3DUpdate<EMFace::Y, S, true>,
  /* Hy */ ReflectingBCUpdate,
  /* Hz */ ReflectingBCUpdate
>;

template<EMSide S>
using PmlXTE = TypeList<
  /* Ex */ ReflectingBCUpdate,
  /* Ey */ Pml3DUpdate<EMFace::X, S, true>,
  /* Ez */ ReflectingBCUpdate,
  /* Hx */ ReflectingBCUpdate,
  /* Hy */ ReflectingBCUpdate,
  /* Hz */ Pml3DUpdate<EMFace::X, S, true>
>;

template<EMSide S>
using PmlYTE = TypeList<
  /* Ex */ Pml3DUpdate<EMFace::X, S, true>,
  /* Ey */ ReflectingBCUpdate,
  /* Ez */ ReflectingBCUpdate,
  /* Hx */ ReflectingBCUpdate,
  /* Hy */ ReflectingBCUpdate,
  /* Hz */ Pml3DUpdate<EMFace::X, S, true>
>;

template<EMSide S>
using PmlX3D = TypeList<
  /* Ex */ ReflectingBCUpdate,
  /* Ey */ Pml3DUpdate<EMFace::X, S, true>,
  /* Ez */ Pml3DUpdate<EMFace::X, S, false>,
  /* Hx */ ReflectingBCUpdate,
  /* Hy */ Pml3DUpdate<EMFace::X, S, false>,
  /* Hz */ Pml3DUpdate<EMFace::X, S, true>
>;

template<EMSide S>
using PmlY3D = TypeList<
  /* Ex */ Pml3DUpdate<EMFace::Y, S, false>,
  /* Ey */ ReflectingBCUpdate,
  /* Ez */ Pml3DUpdate<EMFace::Y, S, true>,
  /* Hx */ Pml3DUpdate<EMFace::Y, S, true>,
  /* Hy */ ReflectingBCUpdate,
  /* Hz */ Pml3DUpdate<EMFace::Y, S, false>
>;

template<EMSide S>
using PmlZ3D = TypeList<
  /* Ex */ Pml3DUpdate<EMFace::Z, S, true>,
  /* Ey */ Pml3DUpdate<EMFace::Z, S, false>,
  /* Ez */ ReflectingBCUpdate,
  /* Hx */ Pml3DUpdate<EMFace::Z, S, false>,
  /* Hy */ Pml3DUpdate<EMFace::Z, S, true>,
  /* Hz */ ReflectingBCUpdate
>;

#endif //BC_DEFINITIONS_H
