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
template<typename T>
using ReflectingData = FaceBCs<
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
using PeriodicData1D = FaceBCs<
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
using PeriodicDataTM = FaceBCs<
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
using PeriodicDataTE = FaceBCs<
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
using PeriodicData3D = FaceBCs<
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
using PmlData1D = FaceBCs<
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
struct PmlDataTMImpl {
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
struct PmlDataTMImpl<T, F, S> {
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
struct PmlDataTMImpl<T, F, S> {
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
using PmlDataTM = typename PmlDataTMImpl<T, F, S>::type;


// -------------------------------------------
// TE PML
template<typename T, EMFace F, EMSide S>
struct PmlDataTEImpl{
  using type = FaceBCs<
    /* Ex */ NullData<EmptyArray2D<T>>,
    /* Ey */ NullData<EmptyArray2D<T>>,
    /* Ez */ NullData<EmptyArray2D<T>>,
    /* Hx */ NullData<EmptyArray2D<T>>,
    /* Hy */ NullData<EmptyArray2D<T>>,
    /* Hz */ NullData<EmptyArray2D<T>>
  >;
};

// TE X-Face
template<typename T, EMFace F, EMSide S>
requires (F == EMFace::X)
struct PmlDataTEImpl<T, F, S> {
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
struct PmlDataTEImpl<T, F, S> {
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
using PmlDataTE = typename PmlDataTEImpl<T, F, S>::type;

// -------------------------------------------
// 3D Pml
template<typename T, EMFace F, EMSide S>
struct PmlData3DImpl;

// 3D X-Face
template<typename T, EMFace F, EMSide S>
requires (F == EMFace::X)
struct PmlData3DImpl<T, F, S> {
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
struct PmlData3DImpl<T, F, S> {
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
struct PmlData3DImpl<T, F, S> {
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
using PmlData3D = typename PmlData3DImpl<T, F, S>::type;

//=================== Boundary Condition Definitions ===================
//======================================================================
using ReflectingBC = TypeList<
  /* Ex */ ReflectingBCUpdate,
  /* Ey */ ReflectingBCUpdate,
  /* Ez */ ReflectingBCUpdate,
  /* Hx */ ReflectingBCUpdate,
  /* Hy */ ReflectingBCUpdate,
  /* Hz */ ReflectingBCUpdate
>;

template<EMFace F, EMSide S>
using Periodic1D = TypeList<
  /* Ex */ ReflectingBCUpdate,
  /* Ey */ ReflectingBCUpdate,
  /* Ez */ Periodic1DUpdate<F, S>,
  /* Hx */ ReflectingBCUpdate,
  /* Hy */ Periodic1DUpdate<F, S>,
  /* Hz */ ReflectingBCUpdate
>;

template<EMFace F, EMSide S>
using PeriodicTM = TypeList<
  /* Ex */ ReflectingBCUpdate,
  /* Ey */ ReflectingBCUpdate,
  /* Ez */ Periodic2DUpdate<F, S>,
  /* Hx */ Periodic2DUpdate<F, S>,
  /* Hy */ Periodic2DUpdate<F, S>,
  /* Hz */ ReflectingBCUpdate
>;

template<EMFace F, EMSide S>
using PeriodicTE = TypeList<
  /* Ex */ Periodic2DUpdate<F, S>,
  /* Ey */ Periodic2DUpdate<F, S>,
  /* Ez */ ReflectingBCUpdate,
  /* Hx */ ReflectingBCUpdate,
  /* Hy */ ReflectingBCUpdate,
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

template<EMFace F, EMSide S>
using Pml1D = TypeList<
  /* Ex */ ReflectingBCUpdate,
  /* Ey */ ReflectingBCUpdate,
  /* Ez */ Pml1DUpdate<F, S>,
  /* Hx */ ReflectingBCUpdate,
  /* Hy */ Pml1DUpdate<F, S>,
  /* Hz */ ReflectingBCUpdate
>;


template<EMFace F, EMSide S>
struct PmlTMImpl {
    using type = TypeList<
    /* Ex */ ReflectingBCUpdate,
    /* Ey */ ReflectingBCUpdate,
    /* Ez */ ReflectingBCUpdate,
    /* Hx */ ReflectingBCUpdate,
    /* Hy */ ReflectingBCUpdate,
    /* Hz */ ReflectingBCUpdate
  >;
};

template<EMFace F, EMSide S>
requires (F == EMFace::X)
struct PmlTMImpl<F, S> {
  using type = TypeList<
    /* Ex */ ReflectingBCUpdate,
    /* Ey */ ReflectingBCUpdate,
    /* Ez */ Pml2DUpdate<EMFace::X, S, false>,
    /* Hx */ ReflectingBCUpdate,
    /* Hy */ Pml2DUpdate<EMFace::X, S, false>,
    /* Hz */ ReflectingBCUpdate
  >;
};

template<EMFace F, EMSide S>
requires (F == EMFace::Y)
struct PmlTMImpl<F, S> {
  using type = TypeList<
    /* Ex */ ReflectingBCUpdate,
    /* Ey */ ReflectingBCUpdate,
    /* Ez */ Pml2DUpdate<EMFace::Y, S, true>,
    /* Hx */ Pml2DUpdate<EMFace::Y, S, true>,
    /* Hy */ ReflectingBCUpdate,
    /* Hz */ ReflectingBCUpdate
  >;
};

// Top-level alias
template<EMFace F, EMSide S>
using PmlTM = typename PmlTMImpl<F, S>::type;

template<EMFace F, EMSide S>
struct PmlTEImpl {
  using type = TypeList<
  /* Ex */ ReflectingBCUpdate,
  /* Ey */ ReflectingBCUpdate,
  /* Ez */ ReflectingBCUpdate,
  /* Hx */ ReflectingBCUpdate,
  /* Hy */ ReflectingBCUpdate,
  /* Hz */ ReflectingBCUpdate
>;
};

template<EMFace F, EMSide S>
requires (F == EMFace::X)
struct PmlTEImpl<F, S> {
  using type = TypeList<
    /* Ex */ ReflectingBCUpdate,
    /* Ey */ Pml2DUpdate<EMFace::X, S, true>,
    /* Ez */ ReflectingBCUpdate,
    /* Hx */ ReflectingBCUpdate,
    /* Hy */ ReflectingBCUpdate,
    /* Hz */ Pml2DUpdate<EMFace::X, S, true>
  >;
};

template<EMFace F, EMSide S>
requires (F == EMFace::Y)
struct PmlTEImpl<F, S> {
  using type = TypeList<
    /* Ex */ Pml2DUpdate<EMFace::Y, S, false>,
    /* Ey */ ReflectingBCUpdate,
    /* Ez */ ReflectingBCUpdate,
    /* Hx */ ReflectingBCUpdate,
    /* Hy */ ReflectingBCUpdate,
    /* Hz */ Pml2DUpdate<EMFace::Y, S, false>
  >;
};

// Top-level alias
template<EMFace F, EMSide S>
using PmlTE = typename PmlTEImpl<F, S>::type;


template<EMFace F, EMSide S>
struct Pml3DImpl;

template<EMFace F, EMSide S>
requires (F == EMFace::X)
struct Pml3DImpl<F, S> {
  using type = TypeList<
    /* Ex */ ReflectingBCUpdate,
    /* Ey */ Pml3DUpdate<EMFace::X, S, true>,
    /* Ez */ Pml3DUpdate<EMFace::X, S, false>,
    /* Hx */ ReflectingBCUpdate,
    /* Hy */ Pml3DUpdate<EMFace::X, S, false>,
    /* Hz */ Pml3DUpdate<EMFace::X, S, true>
  >;
};

template<EMFace F, EMSide S>
requires (F == EMFace::Y)
struct Pml3DImpl<F, S> {
  using type = TypeList<
    /* Ex */ Pml3DUpdate<EMFace::Y, S, false>,
    /* Ey */ ReflectingBCUpdate,
    /* Ez */ Pml3DUpdate<EMFace::Y, S, true>,
    /* Hx */ Pml3DUpdate<EMFace::Y, S, true>,
    /* Hy */ ReflectingBCUpdate,
    /* Hz */ Pml3DUpdate<EMFace::Y, S, false>
  >;
};

template<EMFace F, EMSide S>
requires (F == EMFace::Z)
struct Pml3DImpl<F, S> {
  using type = TypeList<
    /* Ex */ Pml3DUpdate<EMFace::Z, S, true>,
    /* Ey */ Pml3DUpdate<EMFace::Z, S, false>,
    /* Ez */ ReflectingBCUpdate,
    /* Hx */ Pml3DUpdate<EMFace::Z, S, false>,
    /* Hy */ Pml3DUpdate<EMFace::Z, S, true>,
    /* Hz */ ReflectingBCUpdate
  >;
};

// Top-level alias
template<EMFace F, EMSide S>
using Pml3D = typename Pml3DImpl<F, S>::type;

#endif //BC_DEFINITIONS_H
