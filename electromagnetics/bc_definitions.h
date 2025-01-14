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
    /* Ey */ NullData<EmptyArray3D<T>>,
    /* Ez */ NullData<EmptyArray3D<T>>,
    /* Hx */ PeriodicData<Array3D<T>, F, S>,
    /* Hy */ PeriodicData<Array3D<T>, F, S>,
    /* Hz */ NullData<EmptyArray3D<T>>
  >;
};

template<typename T, EMFace F, EMSide S>
requires (F == EMFace::Y)
struct PeriodicData3DImpl<T, F, S> {
  using type = FaceBCs<
    /* Ex */ NullData<EmptyArray3D<T>>,
    /* Ey */ NullData<EmptyArray3D<T>>,
    /* Ez */ NullData<EmptyArray3D<T>>,
    /* Hx */ PeriodicData<Array3D<T>, F, S>,
    /* Hy */ NullData<EmptyArray3D<T>>,
    /* Hz */ PeriodicData<Array3D<T>, F, S>
  >;
};

template<typename T, EMFace F, EMSide S>
requires (F == EMFace::Z)
struct PeriodicData3DImpl<T, F, S> {
  using type = FaceBCs<
    /* Ex */ NullData<EmptyArray3D<T>>,
    /* Ey */ NullData<EmptyArray3D<T>>,
    /* Ez */ NullData<EmptyArray3D<T>>,
    /* Hx */ PeriodicData<Array3D<T>, F, S>,
    /* Hy */ PeriodicData<Array3D<T>, F, S>,
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
    /* Ey */ PMLData<Array3D<T>, F, S, true>,
    /* Ez */ PMLData<Array3D<T>, F, S, true>,
    /* Hx */ NullData<EmptyArray3D<T>>,
    /* Hy */ PMLData<Array3D<T>, F, S, false>,
    /* Hz */ PMLData<Array3D<T>, F, S, false>
  >;
};

// 3D Y-Face
template<typename T, EMFace F, EMSide S>
requires (F == EMFace::Y)
struct PmlData3DImpl<T, F, S> {
  using type = FaceBCs<
    /* Ex */ PMLData<Array3D<T>, F, S, true>,
    /* Ey */ NullData<EmptyArray3D<T>>,
    /* Ez */ PMLData<Array3D<T>, F, S, true>,
    /* Hx */ PMLData<Array3D<T>, F, S, false>,
    /* Hy */ NullData<EmptyArray3D<T>>,
    /* Hz */ PMLData<Array3D<T>, F, S, false>
  >;
};

// 3D Z-Face
template<typename T, EMFace F, EMSide S>
requires (F == EMFace::Z)
struct PmlData3DImpl<T, F, S> {
  using type = FaceBCs<
    /* Ex */ PMLData<Array3D<T>, F, S, true>,
    /* Ey */ PMLData<Array3D<T>, F, S, true>,
    /* Ez */ NullData<EmptyArray3D<T>>,
    /* Hx */ PMLData<Array3D<T>, F, S, false>,
    /* Hy */ PMLData<Array3D<T>, F, S, false>,
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
    /* Ey */ ReflectingBCUpdate,
    /* Ez */ ReflectingBCUpdate,
    /* Hx */ Periodic3DUpdate<F, S>,
    /* Hy */ Periodic3DUpdate<F, S>,
    /* Hz */ ReflectingBCUpdate
  >;
};

template<EMFace F, EMSide S>
requires (F == EMFace::Y)
struct Periodic3DImpl<F, S> {
  using type = TypeList<
    /* Ex */ ReflectingBCUpdate,
    /* Ey */ ReflectingBCUpdate,
    /* Ez */ ReflectingBCUpdate,
    /* Hx */ Periodic3DUpdate<F, S>,
    /* Hy */ ReflectingBCUpdate,
    /* Hz */ Periodic3DUpdate<F, S>
  >;
};

template<EMFace F, EMSide S>
requires (F == EMFace::Z)
struct Periodic3DImpl<F, S> {
  using type = TypeList<
    /* Ex */ ReflectingBCUpdate,
    /* Ey */ ReflectingBCUpdate,
    /* Ez */ ReflectingBCUpdate,
    /* Hx */ Periodic3DUpdate<F, S>,
    /* Hy */ Periodic3DUpdate<F, S>,
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
    /* Ey */ Pml3DUpdate<EMFace::X, S, Derivative::DX, true, true>,
    /* Ez */ Pml3DUpdate<EMFace::X, S, Derivative::DX, false, true>,
    /* Hx */ ReflectingBCUpdate,
    /* Hy */ Pml3DUpdate<EMFace::X, S, Derivative::DX, false, false>,
    /* Hz */ Pml3DUpdate<EMFace::X, S, Derivative::DX,  true, false>
  >;
};

template<EMFace F, EMSide S>
requires (F == EMFace::Y)
struct Pml3DImpl<F, S> {
  using type = TypeList<
    /* Ex */ Pml3DUpdate<EMFace::Y, S, Derivative::DY, false, true>,
    /* Ey */ ReflectingBCUpdate,
    /* Ez */ Pml3DUpdate<EMFace::Y, S, Derivative::DY, true, true>,
    /* Hx */ Pml3DUpdate<EMFace::Y, S, Derivative::DY, true, false>,
    /* Hy */ ReflectingBCUpdate,
    /* Hz */ Pml3DUpdate<EMFace::Y, S, Derivative::DY, false, false>
  >;
};

template<EMFace F, EMSide S>
requires (F == EMFace::Z)
struct Pml3DImpl<F, S> {
  using type = TypeList<
    /* Ex */ Pml3DUpdate<EMFace::Z, S, Derivative::DZ, true, true>,
    /* Ey */ Pml3DUpdate<EMFace::Z, S, Derivative::DZ, false, true>,
    /* Ez */ ReflectingBCUpdate,
    /* Hx */ Pml3DUpdate<EMFace::Z, S, Derivative::DZ, false, false>,
    /* Hy */ Pml3DUpdate<EMFace::Z, S, Derivative::DZ, true, false>,
    /* Hz */ ReflectingBCUpdate
  >;
};

// Top-level alias
template<EMFace F, EMSide S>
using Pml3D = typename Pml3DImpl<F, S>::type;

#endif // BC_DEFINITIONS_H
