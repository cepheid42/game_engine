//
// Created by cepheid on 10/17/24.
//

#ifndef ELECTROMAGNETICS_H
#define ELECTROMAGNETICS_H

#include "em_definitions.h"
#include "bc_definitions.h"

inline constexpr size_t SELECT_EMDATA = 4; // todo: these can be combined into one value?
inline constexpr size_t SELECT_EMSOLVER = 4;

inline constexpr size_t SELECT_BCDATA[6] = {8, 0, 11, 0, 12, 0}; // Xlo, Xhi, Ylo, Yhi, Zlo, Zhi
inline constexpr size_t SELECT_BCSOLVER[6] = {6, 0, 6, 0, 6, 0}; // Xlo, Xhi, Ylo, Yhi, Zlo, Zhi

//=================== Boundary Condition Selectors ===================
//====================================================================
template<typename T, EMFace F, EMSide S>
using BCDataTL = TypeList<
  Reflecting_Face<T, S>,    // 0
  Periodic_Face1D<T, F, S>, // 1
  Periodic_FaceTM<T, F, S>, // 2
  Periodic_FaceTE<T, F, S>, // 3
  Periodic_Face3D<T, F, S>, // 4
  PML_XFace1D<T, S>,        // 5
  PML_XFaceTM<T, S>,        // 6
  PML_XFaceTE<T, S>,        // 7
  PML_XFace3D<T, S>,        // 8
  PML_YFaceTM<T, S>,        // 9
  PML_YFaceTE<T, S>,        // 10
  PML_YFace3D<T, S>,        // 11
  PML_ZFace3D<T, S>         // 12
>;

template<typename T>
using bcdata_t = BCData<
  TypeListAt<SELECT_BCDATA[0], BCDataTL<T, EMFace::X, EMSide::Lo>>,
  TypeListAt<SELECT_BCDATA[1], BCDataTL<T, EMFace::X, EMSide::Hi>>,
  TypeListAt<SELECT_BCDATA[2], BCDataTL<T, EMFace::Y, EMSide::Lo>>,
  TypeListAt<SELECT_BCDATA[3], BCDataTL<T, EMFace::Y, EMSide::Hi>>,
  TypeListAt<SELECT_BCDATA[4], BCDataTL<T, EMFace::Z, EMSide::Lo>>,
  TypeListAt<SELECT_BCDATA[5], BCDataTL<T, EMFace::Z, EMSide::Hi>>
>;

template<EMFace F, EMSide S, bool Negate>
using BCTypeTL = TypeList<
  ReflectingBC,
  Periodic1D,
  Periodic2D<F, S>,
  Periodic3D<F, S>,
  Pml1D<S>,
  Pml2D<F, S, Negate>,
  Pml3D<F, S, Negate>
>;

template<size_t I, EMFace F, EMSide S, bool Negate>
using BCType = TypeListAt<I, BCTypeTL<F, S, Negate>>;

/*
 *
 * todo: figure out how to group the functions for the faces together
 *
 */
//=================== Electromagnetics Selectors ===================
//==================================================================
template<typename T>
using EMDataTL = TypeList<
  emdataNone<T>, // 0
  emdata1D<T>,   // 1
  emdataTM<T>,   // 2
  emdataTE<T>,   // 3
  emdata3D<T>    // 4
>; // Typelist for choosing type of EMData

template<typename T>
using emdata_t = TypeListAt<SELECT_EMDATA, EMDataTL<T>>;

template<typename T>
using EMTypeTL = TypeList<
  EMNull<T>, // 0
  EM1D<T>,   // 1
  EMTM<T>,   // 2
  EMTE<T>,   // 3
  EM3D<T>    // 4
>; // Typelist (of typelists) for choosing type of EM Solver

template<typename T>
using EMType = TypeListAt<SELECT_EMSOLVER, EMTypeTL<emdata_t<T>>>; // Selects desired typelist of integrators

template<typename T>
using EMSolver = Electromagnetics<
  TypeListAt<0, EMType<T>>, // Ex
  TypeListAt<1, EMType<T>>, // Ey
  TypeListAt<2, EMType<T>>, // Ez
  TypeListAt<3, EMType<T>>, // Hx
  TypeListAt<4, EMType<T>>, // Hy
  TypeListAt<5, EMType<T>>  // Hz
>; // Auto-selects full BC's per face and EMSolver per field component based on chosen settings.

#endif //ELECTROMAGNETICS_H
