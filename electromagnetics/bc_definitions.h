//
// Created by cepheid on 11/11/24.
//

#ifndef BC_DEFINITIONS_H
#define BC_DEFINITIONS_H

#include "bc_data.h"
#include "boundaries.h"

// todo: Does this need to have namespaces? Are these usings polluting upstream?
using tf::electromagnetics::traits::EMFace;
using tf::electromagnetics::traits::EMSide;
using tf::electromagnetics::boundaries::NullData;
using tf::electromagnetics::boundaries::FaceBCs;
using tf::electromagnetics::boundaries::PeriodicData;
using tf::electromagnetics::boundaries::PMLData;
using tf::electromagnetics::boundaries::ReflectingBCUpdate;
using tf::electromagnetics::boundaries::Periodic3DUpdate;
using tf::electromagnetics::boundaries::Pml3DUpdate;
using tf::types::Array3D;
using tf::types::EmptyArray3D;

//=================== BCData Definitions ===================
//==========================================================
// -------------------------------------------
// Reflecting Boundary
template<typename T>
using ReflectingData = FaceBCs<
  /* Ex */ NullData<EmptyArray3D<T>>,
  /* Ey */ NullData<EmptyArray3D<T>>,
  /* Ez */ NullData<EmptyArray3D<T>>,
  /* Hx */ NullData<EmptyArray3D<T>>,
  /* Hy */ NullData<EmptyArray3D<T>>,
  /* Hz */ NullData<EmptyArray3D<T>>
>;

// -------------------------------------------
// 3D Periodic Boundary
template<typename T, EMFace F, EMSide S>
struct PeriodicData3DImpl;

template<typename T, EMFace F, EMSide S>
requires (F == EMFace::X)
struct PeriodicData3DImpl<T, F, S> {
  using type = FaceBCs<
    /* Ex */ NullData<EmptyArray3D<T>>,
    /* Ey */ PeriodicData<Array3D<T>, F, S>,
    /* Ez */ PeriodicData<Array3D<T>, F, S>,
    /* Hx */ NullData<EmptyArray3D<T>>,
    /* Hy */ NullData<EmptyArray3D<T>>,
    /* Hz */ NullData<EmptyArray3D<T>>
  >;
};

template<typename T, EMFace F, EMSide S>
requires (F == EMFace::Y)
struct PeriodicData3DImpl<T, F, S> {
  using type = FaceBCs<
    /* Ex */ PeriodicData<Array3D<T>, F, S>,
    /* Ey */ NullData<EmptyArray3D<T>>,
    /* Ez */ PeriodicData<Array3D<T>, F, S>,
    /* Hx */ NullData<EmptyArray3D<T>>,
    /* Hy */ NullData<EmptyArray3D<T>>,
    /* Hz */ NullData<EmptyArray3D<T>>
  >;
};

template<typename T, EMFace F, EMSide S>
requires (F == EMFace::Z)
struct PeriodicData3DImpl<T, F, S> {
  using type = FaceBCs<
    /* Ex */ PeriodicData<Array3D<T>, F, S>,
    /* Ey */ PeriodicData<Array3D<T>, F, S>,
    /* Ez */ NullData<EmptyArray3D<T>>,
    /* Hx */ NullData<EmptyArray3D<T>>,
    /* Hy */ NullData<EmptyArray3D<T>>,
    /* Hz */ NullData<EmptyArray3D<T>>
  >;
};

// Top-level alias
template<typename T, EMFace F, EMSide S>
using PeriodicData3D = typename PeriodicData3DImpl<T, F, S>::type;


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

// ===============================================
template<EMFace F, EMSide S>
struct Periodic3DImpl;

template<EMFace F, EMSide S>
requires (F == EMFace::X)
struct Periodic3DImpl<F, S> {
  using type = TypeList<
    /* Ex */ ReflectingBCUpdate,
    /* Ey */ Periodic3DUpdate<F, S>,
    /* Ez */ Periodic3DUpdate<F, S>,
    /* Hx */ ReflectingBCUpdate,
    /* Hy */ ReflectingBCUpdate,
    /* Hz */ ReflectingBCUpdate
  >;
};

template<EMFace F, EMSide S>
requires (F == EMFace::Y)
struct Periodic3DImpl<F, S> {
  using type = TypeList<
    /* Ex */ Periodic3DUpdate<F, S>,
    /* Ey */ ReflectingBCUpdate,
    /* Ez */ Periodic3DUpdate<F, S>,
    /* Hx */ ReflectingBCUpdate,
    /* Hy */ ReflectingBCUpdate,
    /* Hz */ ReflectingBCUpdate
  >;
};

template<EMFace F, EMSide S>
requires (F == EMFace::Z)
struct Periodic3DImpl<F, S> {
  using type = TypeList<
    /* Ex */ Periodic3DUpdate<F, S>,
    /* Ey */ Periodic3DUpdate<F, S>,
    /* Ez */ ReflectingBCUpdate,
    /* Hx */ ReflectingBCUpdate,
    /* Hy */ ReflectingBCUpdate,
    /* Hz */ ReflectingBCUpdate
  >;
};

// Top-level alias
template<EMFace F, EMSide S>
using Periodic3D = typename Periodic3DImpl<F, S>::type;

// ===============================================
template<EMFace F, EMSide S>
struct Pml3DImpl;

template<EMFace F, EMSide S>
requires (F == EMFace::X)
struct Pml3DImpl<F, S> {
  using type = TypeList<
    /* Ex */ ReflectingBCUpdate,
    /* Ey */ Pml3DUpdate<EMFace::X, S, Derivative::DX, true, false>,
    /* Ez */ Pml3DUpdate<EMFace::X, S, Derivative::DX, false, false>,
    /* Hx */ ReflectingBCUpdate,
    /* Hy */ Pml3DUpdate<EMFace::X, S, Derivative::DX, false, true>,
    /* Hz */ Pml3DUpdate<EMFace::X, S, Derivative::DX,  true, true>
  >;
};

template<EMFace F, EMSide S>
requires (F == EMFace::Y)
struct Pml3DImpl<F, S> {
  using type = TypeList<
    /* Ex */ Pml3DUpdate<EMFace::Y, S, Derivative::DY, false, false>,
    /* Ey */ ReflectingBCUpdate,
    /* Ez */ Pml3DUpdate<EMFace::Y, S, Derivative::DY, true, false>,
    /* Hx */ Pml3DUpdate<EMFace::Y, S, Derivative::DY, true, true>,
    /* Hy */ ReflectingBCUpdate,
    /* Hz */ Pml3DUpdate<EMFace::Y, S, Derivative::DY, false, true>
  >;
};

template<EMFace F, EMSide S>
requires (F == EMFace::Z)
struct Pml3DImpl<F, S> {
  using type = TypeList<
    /* Ex */ Pml3DUpdate<EMFace::Z, S, Derivative::DZ, true, false>,
    /* Ey */ Pml3DUpdate<EMFace::Z, S, Derivative::DZ, false, false>,
    /* Ez */ ReflectingBCUpdate,
    /* Hx */ Pml3DUpdate<EMFace::Z, S, Derivative::DZ, false, true>,
    /* Hy */ Pml3DUpdate<EMFace::Z, S, Derivative::DZ, true, true>,
    /* Hz */ ReflectingBCUpdate
  >;
};

// Top-level alias
template<EMFace F, EMSide S>
using Pml3D = typename Pml3DImpl<F, S>::type;

#endif // BC_DEFINITIONS_H
