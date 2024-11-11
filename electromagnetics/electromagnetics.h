//
// Created by cepheid on 10/17/24.
//

#ifndef ELECTROMAGNETICS_H
#define ELECTROMAGNETICS_H

#include "em_definitions.h"
#include "bc_definitions.h"

inline constexpr size_t SELECT_EMDATA = 1; // todo: these can be combined into one value?
inline constexpr size_t SELECT_EMSOLVER = 1;

static constexpr size_t SELECT_BCDATA[6] = {8, 12, 12, 12, 12, 12}; // Xlo, Xhi, Ylo, Yhi, Zlo, Zhi
static constexpr size_t SELECT_BCSOLVER[6] = {1, 0, 0, 0, 0, 0}; // Xlo, Xhi, Ylo, Yhi, Zlo, Zhi

// using fp_t = double;
// constexpr size_t DIM = 1;
// constexpr fp_t cfl = 0.95 / std::sqrt(static_cast<fp_t>(DIM));

//=================== Boundary Condition Selectors ===================
//====================================================================
template<typename T, EMFace F, EMSide S>
using BCFaceTL = TypeList<
  PML_XFace1D<T, S>,        // 0
  PML_XFaceTM<T, S>,        // 1
  PML_XFaceTE<T, S>,        // 2
  PML_XFace3D<T, S>,        // 3
  PML_YFaceTM<T, S>,        // 4
  PML_YFaceTE<T, S>,        // 5
  PML_YFace3D<T, S>,        // 6
  PML_ZFace3D<T, S>,        // 7
  Periodic_Face1D<T, F, S>, // 8
  Periodic_FaceTM<T, F, S>, // 9
  Periodic_FaceTE<T, F, S>, // 10
  Periodic_Face3D<T, F, S>, // 11
  Reflecting_Face<T, S>     // 12
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
  BCReflecting<T>,
  BC1DPeriodic<T>,
  BC1DPml<T>
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
  TypeListAt<4, BCType<SELECT_BCSOLVER[0], T>> // X0 Face
>; // Auto-selects full BC's per face and EMSolver per field component based on chosen settings.

#endif //ELECTROMAGNETICS_H
