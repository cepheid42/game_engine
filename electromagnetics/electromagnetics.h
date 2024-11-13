//
// Created by cepheid on 10/17/24.
//

#ifndef ELECTROMAGNETICS_H
#define ELECTROMAGNETICS_H

#include "em_definitions.h"
#include "bc_definitions.h"

inline constexpr size_t SELECT_EMDATA = 2; // todo: these can be combined into one value?
inline constexpr size_t SELECT_EMSOLVER = 2;

static constexpr size_t SELECT_BCDATA[6] = {2, 0, 0, 0, 0, 0}; // Xlo, Xhi, Ylo, Yhi, Zlo, Zhi
static constexpr size_t SELECT_BCSOLVER[6] = {2, 0, 2, 0, 0, 0}; // Xlo, Xhi, Ylo, Yhi, Zlo, Zhi

//=================== Boundary Condition Selectors ===================
//====================================================================
template<typename T, EMFace F, EMSide S>
using BCFaceTL = TypeList<
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
  TypeListAt<SELECT_BCDATA[0], BCFaceTL<T, EMFace::X, EMSide::Lo>>,
  TypeListAt<SELECT_BCDATA[1], BCFaceTL<T, EMFace::X, EMSide::Hi>>,
  TypeListAt<SELECT_BCDATA[2], BCFaceTL<T, EMFace::Y, EMSide::Lo>>,
  TypeListAt<SELECT_BCDATA[3], BCFaceTL<T, EMFace::Y, EMSide::Hi>>,
  TypeListAt<SELECT_BCDATA[4], BCFaceTL<T, EMFace::Z, EMSide::Lo>>,
  TypeListAt<SELECT_BCDATA[5], BCFaceTL<T, EMFace::Z, EMSide::Hi>>
>;


template<typename T>
using BCTypeTL = TypeList<
  BCReflecting<T>,  // 0
  BCPeriodic1D<T>,  // 1
  BCPeriodicTM<T>,  // 2
  BCPeriodicTE<T>,  // 3
  BCPeriodic3D<T>,  // 4
  BCPml1D<T>,       // 5
  BCPmlTM_XFace<T, EMFace::X>, // 6
  BCPmlTM_YFace<T, EMFace::Y>, // 7
  BCPmlTE_XFace<T, EMFace::X>, // 8
  BCPmlTE_YFace<T, EMFace::Y>, // 9
  BCPml3D_XFace<T, EMFace::X>, // 10
  BCPml3D_YFace<T, EMFace::Y>, // 11
  BCPml3D_ZFace<T, EMFace::Z>  // 12
>;

template<size_t I, typename T>
using BCType = TypeListAt<I, BCTypeTL<bcdata_t<T>>>;


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
  TypeListAt<5, EMType<T>>, // Hz,
  BCType<SELECT_BCSOLVER[0], T>, // X0 Face
  BCType<SELECT_BCSOLVER[1], T>, // X1 Face
  BCType<SELECT_BCSOLVER[2], T>, // Y0 Face
  BCType<SELECT_BCSOLVER[3], T>, // Y1 Face
  BCType<SELECT_BCSOLVER[4], T>, // Z0 Face
  BCType<SELECT_BCSOLVER[5], T>  // Z1 Face
>; // Auto-selects full BC's per face and EMSolver per field component based on chosen settings.

#endif //ELECTROMAGNETICS_H
