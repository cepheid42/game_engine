//
// Created by cepheid on 10/17/24.
//

#ifndef ELECTROMAGNETICS_H
#define ELECTROMAGNETICS_H

// #include "electromagnetics.param"
#include "em_definitions.h"
#include "bc_definitions.h"

// todo: Does this need to have namespaces?
using tf::electromagnetics::EMData;
using tf::electromagnetics::Electromagnetics;
using tf::electromagnetics::boundaries::BCData;
using tf::electromagnetics::boundaries::Boundary;


// todo: This checks for invalid BC combos that will compile & run but eventually blow up.
//        Probably can be made more comprehensive later with concepts or something.
static constexpr bool valid_x_combo = SELECT_BCSOLVER[0] + 4 != SELECT_BCSOLVER[1];
static constexpr bool valid_y_combo = SELECT_BCSOLVER[2] + 4 != SELECT_BCSOLVER[3];
static constexpr bool valid_z_combo = SELECT_BCSOLVER[4] + 4 != SELECT_BCSOLVER[5];
static_assert(valid_x_combo and valid_y_combo and valid_z_combo,
  "Periodic lower boundary and PML upper boundary is an invalid combination.");

//=================== Boundary Condition Selectors ===================
//====================================================================
template<typename T, EMFace F, EMSide S>
using BCDataTL = TypeList<
  ReflectingData<T>,        // 0
  PeriodicData1D<T, F, S>,  // 1
  PeriodicDataTM<T, F, S>,  // 2
  PeriodicDataTE<T, F, S>,  // 3
  PeriodicData3D<T, F, S>,  // 4
  PmlData1D<T, F, S>,       // 5
  PmlDataTM<T, F, S>,       // 6
  PmlDataTE<T, F, S>,       // 7
  PmlData3D<T, F, S>,       // 8
  PmlData2D6C<T, F, S>      // 9
>;

template<typename T>
using bcdata_t = BCData<
  TypeListAt<SELECT_BCSOLVER[0], BCDataTL<T, EMFace::X, EMSide::Lo>>,
  TypeListAt<SELECT_BCSOLVER[1], BCDataTL<T, EMFace::X, EMSide::Hi>>,
  TypeListAt<SELECT_BCSOLVER[2], BCDataTL<T, EMFace::Y, EMSide::Lo>>,
  TypeListAt<SELECT_BCSOLVER[3], BCDataTL<T, EMFace::Y, EMSide::Hi>>,
  TypeListAt<SELECT_BCSOLVER[4], BCDataTL<T, EMFace::Z, EMSide::Lo>>,
  TypeListAt<SELECT_BCSOLVER[5], BCDataTL<T, EMFace::Z, EMSide::Hi>>
>;

template<EMFace F, EMSide S>
using BCTypeTL = TypeList<
  ReflectingBC,     // 0
  Periodic1D<F, S>, // 1
  PeriodicTM<F, S>, // 2
  PeriodicTE<F, S>, // 3
  Periodic3D<F, S>, // 4
  Pml1D<F, S>,      // 5
  PmlTM<F, S>,      // 6
  PmlTE<F, S>,      // 7
  Pml3D<F, S>,      // 8
  Pml2D6C<F, S>     // 9
>;

template<size_t I, EMFace F, EMSide S>
using BCType = TypeListAt<I, BCTypeTL<F, S>>;

template<typename TL>
struct boundary_t {
  using Ex = Boundary<TypeListAt<0, TL>>;
  using Ey = Boundary<TypeListAt<1, TL>>;
  using Ez = Boundary<TypeListAt<2, TL>>;
  using Hx = Boundary<TypeListAt<3, TL>>;
  using Hy = Boundary<TypeListAt<4, TL>>;
  using Hz = Boundary<TypeListAt<5, TL>>;
};

using bcx0_t = boundary_t<BCType<SELECT_BCSOLVER[0], EMFace::X, EMSide::Lo>>;
using bcx1_t = boundary_t<BCType<SELECT_BCSOLVER[1], EMFace::X, EMSide::Hi>>;
using bcy0_t = boundary_t<BCType<SELECT_BCSOLVER[2], EMFace::Y, EMSide::Lo>>;
using bcy1_t = boundary_t<BCType<SELECT_BCSOLVER[3], EMFace::Y, EMSide::Hi>>;
using bcz0_t = boundary_t<BCType<SELECT_BCSOLVER[4], EMFace::Z, EMSide::Lo>>;
using bcz1_t = boundary_t<BCType<SELECT_BCSOLVER[5], EMFace::Z, EMSide::Hi>>;

//=================== Electromagnetics Selectors ===================
//==================================================================
template<typename T>
using EMDataTL = TypeList<
  emdataNone<T>, // 0
  emdata1D<T>,   // 1
  emdataTMz<T>,   // 2
  emdataTEy<T>,   // 3
  emdata3D<T>,   // 4
  emdata2D6C<T>  // 5
>; // Typelist for choosing type of EMData

template<typename T>
using emdata_t = TypeListAt<SELECT_EMSOLVER, EMDataTL<T>>;

template<typename T>
using EMTypeTL = TypeList<
  EMNull<T>, // 0
  EM1D<T>,   // 1
  EMTMz<T>,   // 2
  EMTEy<T>,   // 3
  EM3D<T>,   // 4
  EM2D6C<T>  // 5
>; // Typelist (of typelists) for choosing type of EM Solver

template<typename T>
using EMType = TypeListAt<SELECT_EMSOLVER, EMTypeTL<emdata_t<T>>>; // Selects desired typelist of integrators

template<typename T> using EXI = TypeListAt<0, EMType<T>>;
template<typename T> using EYI = TypeListAt<1, EMType<T>>;
template<typename T> using EZI = TypeListAt<2, EMType<T>>;
template<typename T> using HXI = TypeListAt<3, EMType<T>>;
template<typename T> using HYI = TypeListAt<4, EMType<T>>;
template<typename T> using HZI = TypeListAt<5, EMType<T>>;

template<typename T>
using EMSolver = Electromagnetics<
  EXI<T>, EYI<T>, EZI<T>,
  HXI<T>, HYI<T>, HZI<T>,
  bcx0_t, bcx1_t,
  bcy0_t, bcy1_t,
  bcz0_t, bcz1_t
>; // Auto-selects full BC's per face and EMSolver per field component based on chosen settings.

#endif //ELECTROMAGNETICS_H
