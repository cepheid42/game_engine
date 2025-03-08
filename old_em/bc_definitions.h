//
// Created by cepheid on 11/11/24.
//

#ifndef BC_DEFINITIONS_H
#define BC_DEFINITIONS_H

#include "bc_data.h"
#include "boundaries.h"

namespace tf::electromagnetics
{
  //=================== BCData Definitions ===================
  //==========================================================
  // -------------------------------------------
  // Reflecting Boundary
  template<typename T>
  using ReflectingData = boundaries::FaceBCs<
    /* Ex */ boundaries::NullData<EmptyArray3D<T>>,
    /* Ey */ boundaries::NullData<EmptyArray3D<T>>,
    /* Ez */ boundaries::NullData<EmptyArray3D<T>>,
    /* Hx */ boundaries::NullData<EmptyArray3D<T>>,
    /* Hy */ boundaries::NullData<EmptyArray3D<T>>,
    /* Hz */ boundaries::NullData<EmptyArray3D<T>>
  >;

  // -------------------------------------------
  // 3D Periodic Boundary
  template<typename T, EMFace F, EMSide S>
  struct PeriodicData3DImpl;

  template<typename T, EMFace F, EMSide S>
  requires (F == EMFace::X)
  struct PeriodicData3DImpl<T, F, S> {
    using type = boundaries::FaceBCs<
      /* Ex */ boundaries::NullData<EmptyArray3D<T>>,
      /* Ey */ boundaries::NullData<EmptyArray3D<T>>,
      /* Ez */ boundaries::NullData<EmptyArray3D<T>>,
      /* Hx */ boundaries::PeriodicData<Array3D<T>, F, S>,
      /* Hy */ boundaries::PeriodicData<Array3D<T>, F, S>,
      /* Hz */ boundaries::NullData<EmptyArray3D<T>>
    >;
  };

  template<typename T, EMFace F, EMSide S>
  requires (F == EMFace::Y)
  struct PeriodicData3DImpl<T, F, S> {
    using type = boundaries::FaceBCs<
      /* Ex */ boundaries::NullData<EmptyArray3D<T>>,
      /* Ey */ boundaries::NullData<EmptyArray3D<T>>,
      /* Ez */ boundaries::NullData<EmptyArray3D<T>>,
      /* Hx */ boundaries::PeriodicData<Array3D<T>, F, S>,
      /* Hy */ boundaries::NullData<EmptyArray3D<T>>,
      /* Hz */ boundaries::PeriodicData<Array3D<T>, F, S>
    >;
  };

  template<typename T, EMFace F, EMSide S>
  requires (F == EMFace::Z)
  struct PeriodicData3DImpl<T, F, S> {
    using type = boundaries::FaceBCs<
      /* Ex */ boundaries::NullData<EmptyArray3D<T>>,
      /* Ey */ boundaries::NullData<EmptyArray3D<T>>,
      /* Ez */ boundaries::NullData<EmptyArray3D<T>>,
      /* Hx */ boundaries::PeriodicData<Array3D<T>, F, S>,
      /* Hy */ boundaries::PeriodicData<Array3D<T>, F, S>,
      /* Hz */ boundaries::NullData<EmptyArray3D<T>>
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
    using type = boundaries::FaceBCs<
      /* Ex */ boundaries::NullData<EmptyArray3D<T>>,
      /* Ey */ boundaries::PMLData<Array3D<T>, F, S, true>,
      /* Ez */ boundaries::PMLData<Array3D<T>, F, S, true>,
      /* Hx */ boundaries::NullData<EmptyArray3D<T>>,
      /* Hy */ boundaries::PMLData<Array3D<T>, F, S, false>,
      /* Hz */ boundaries::PMLData<Array3D<T>, F, S, false>
    >;
  };

  // 3D Y-Face
  template<typename T, EMFace F, EMSide S>
  requires (F == EMFace::Y)
  struct PmlData3DImpl<T, F, S> {
    using type = boundaries::FaceBCs<
      /* Ex */ boundaries::PMLData<Array3D<T>, F, S, true>,
      /* Ey */ boundaries::NullData<EmptyArray3D<T>>,
      /* Ez */ boundaries::PMLData<Array3D<T>, F, S, true>,
      /* Hx */ boundaries::PMLData<Array3D<T>, F, S, false>,
      /* Hy */ boundaries::NullData<EmptyArray3D<T>>,
      /* Hz */ boundaries::PMLData<Array3D<T>, F, S, false>
    >;
  };

  // 3D Z-Face
  template<typename T, EMFace F, EMSide S>
  requires (F == EMFace::Z)
  struct PmlData3DImpl<T, F, S> {
    using type = boundaries::FaceBCs<
      /* Ex */ boundaries::PMLData<Array3D<T>, F, S, true>,
      /* Ey */ boundaries::PMLData<Array3D<T>, F, S, true>,
      /* Ez */ boundaries::NullData<EmptyArray3D<T>>,
      /* Hx */ boundaries::PMLData<Array3D<T>, F, S, false>,
      /* Hy */ boundaries::PMLData<Array3D<T>, F, S, false>,
      /* Hz */ boundaries::NullData<EmptyArray3D<T>>
    >;
  };

  // Top-level alias
  template<typename T, EMFace F, EMSide S>
  using PmlData3D = typename PmlData3DImpl<T, F, S>::type;


  //=================== Boundary Condition Definitions ===================
  //======================================================================
  using ReflectingBC = TypeList<
    /* Ex */ boundaries::ReflectingBCUpdate,
    /* Ey */ boundaries::ReflectingBCUpdate,
    /* Ez */ boundaries::ReflectingBCUpdate,
    /* Hx */ boundaries::ReflectingBCUpdate,
    /* Hy */ boundaries::ReflectingBCUpdate,
    /* Hz */ boundaries::ReflectingBCUpdate
  >;

  // ===============================================
  template<EMFace F, EMSide S>
  struct Periodic3DImpl;

  template<EMFace F, EMSide S>
  requires (F == EMFace::X)
  struct Periodic3DImpl<F, S> {
    using type = TypeList<
      /* Ex */ boundaries::ReflectingBCUpdate,
      /* Ey */ boundaries::ReflectingBCUpdate,
      /* Ez */ boundaries::ReflectingBCUpdate,
      /* Hx */ boundaries::Periodic3DUpdate<F, S>,
      /* Hy */ boundaries::Periodic3DUpdate<F, S>,
      /* Hz */ boundaries::ReflectingBCUpdate
    >;
  };

  template<EMFace F, EMSide S>
  requires (F == EMFace::Y)
  struct Periodic3DImpl<F, S> {
    using type = TypeList<
      /* Ex */ boundaries::ReflectingBCUpdate,
      /* Ey */ boundaries::ReflectingBCUpdate,
      /* Ez */ boundaries::ReflectingBCUpdate,
      /* Hx */ boundaries::Periodic3DUpdate<F, S>,
      /* Hy */ boundaries::ReflectingBCUpdate,
      /* Hz */ boundaries::Periodic3DUpdate<F, S>
    >;
  };

  template<EMFace F, EMSide S>
  requires (F == EMFace::Z)
  struct Periodic3DImpl<F, S> {
    using type = TypeList<
      /* Ex */ boundaries::ReflectingBCUpdate,
      /* Ey */ boundaries::ReflectingBCUpdate,
      /* Ez */ boundaries::ReflectingBCUpdate,
      /* Hx */ boundaries::Periodic3DUpdate<F, S>,
      /* Hy */ boundaries::Periodic3DUpdate<F, S>,
      /* Hz */ boundaries::ReflectingBCUpdate
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
      /* Ex */ boundaries::ReflectingBCUpdate,
      /* Ey */ boundaries::Pml3DUpdate<EMFace::X, S, Derivative::DX, true, true>,
      /* Ez */ boundaries::Pml3DUpdate<EMFace::X, S, Derivative::DX, false, true>,
      /* Hx */ boundaries::ReflectingBCUpdate,
      /* Hy */ boundaries::Pml3DUpdate<EMFace::X, S, Derivative::DX, false, false>,
      /* Hz */ boundaries::Pml3DUpdate<EMFace::X, S, Derivative::DX,  true, false>
    >;
  };

  template<EMFace F, EMSide S>
  requires (F == EMFace::Y)
  struct Pml3DImpl<F, S> {
    using type = TypeList<
      /* Ex */ boundaries::Pml3DUpdate<EMFace::Y, S, Derivative::DY, false, true>,
      /* Ey */ boundaries::ReflectingBCUpdate,
      /* Ez */ boundaries::Pml3DUpdate<EMFace::Y, S, Derivative::DY, true, true>,
      /* Hx */ boundaries::Pml3DUpdate<EMFace::Y, S, Derivative::DY, true, false>,
      /* Hy */ boundaries::ReflectingBCUpdate,
      /* Hz */ boundaries::Pml3DUpdate<EMFace::Y, S, Derivative::DY, false, false>
    >;
  };

  template<EMFace F, EMSide S>
  requires (F == EMFace::Z)
  struct Pml3DImpl<F, S> {
    using type = TypeList<
      /* Ex */ boundaries::Pml3DUpdate<EMFace::Z, S, Derivative::DZ, true, true>,
      /* Ey */ boundaries::Pml3DUpdate<EMFace::Z, S, Derivative::DZ, false, true>,
      /* Ez */ boundaries::ReflectingBCUpdate,
      /* Hx */ boundaries::Pml3DUpdate<EMFace::Z, S, Derivative::DZ, false, false>,
      /* Hy */ boundaries::Pml3DUpdate<EMFace::Z, S, Derivative::DZ, true, false>,
      /* Hz */ boundaries::ReflectingBCUpdate
    >;
  };

  // Top-level alias
  template<EMFace F, EMSide S>
  using Pml3D = typename Pml3DImpl<F, S>::type;

} // end namespace tf::old_em

#endif // BC_DEFINITIONS_H
