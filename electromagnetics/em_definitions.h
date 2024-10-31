//
// Created by cepheid on 10/17/24.
//

#ifndef EM_DEFINITIONS_H
#define EM_DEFINITIONS_H

#include "core/typelist.h"
#include "em_data.h"
#include "em_solver.h"
#include "bc_data.h"
#include "boundaries.h"
#include "offsets.h"

//=================== EMData Definitions ========================
//===============================================================
template<typename Array>
struct enabled {
  using value_t = typename Array::value_t;
  using dimension_t = typename Array::dimension_t;
  using array_t = Array;
};

template<typename Array>
struct disabled {
  using value_t = typename Array::value_t;
  using dimension_t = typename Array::dimension_t;
  using array_t = Array;
};

template<typename T>
using emdataNone = EMData<
  disabled<EmptyArray1D<T>>, // Ex
  disabled<EmptyArray1D<T>>, // Ey
  disabled<EmptyArray1D<T>>, // Ez
  disabled<EmptyArray1D<T>>, // Hx
  disabled<EmptyArray1D<T>>, // Hy
  disabled<EmptyArray1D<T>>  // Hz
>;

template<typename T>
using emdata1D = EMData<
  disabled<EmptyArray1D<T>>, // Ex
  disabled<EmptyArray1D<T>>, // Ey
  enabled<Array1D<T>>,       // Ez
  disabled<EmptyArray1D<T>>, // Hx
  enabled<Array1D<T>>,       // Hy
  disabled<EmptyArray1D<T>>  // Hz
>;

template<typename T>
using emdataTM = EMData<
  disabled<EmptyArray2D<T>>, // Ex
  disabled<EmptyArray2D<T>>, // Ey
  enabled<Array2D<T>>,       // Ez
  enabled<Array2D<T>>,       // Hx
  enabled<Array2D<T>>,       // Hy
  disabled<EmptyArray2D<T>>  // Hz
>;

template<typename T>
using emdataTE = EMData<
  enabled<Array2D<T>>,       // Ex
  enabled<Array2D<T>>,       // Ey
  disabled<EmptyArray2D<T>>, // Ez
  disabled<EmptyArray2D<T>>, // Hx
  disabled<EmptyArray2D<T>>, // Hy
  enabled<Array2D<T>>        // Hz
>;

template<typename T>
using emdata3D = EMData<
  enabled<Array3D<T>>, // Ex
  enabled<Array3D<T>>, // Ey
  enabled<Array3D<T>>, // Ez
  enabled<Array3D<T>>, // Hx
  enabled<Array3D<T>>, // Hy
  enabled<Array3D<T>>  // Hz
>;

//=================== Electromagnetics Definitions ========================
//=========================================================================
template<typename T>
using EMNull = TypeList<
  /* Ex */ FieldIntegratorNull<T>,
  /* Ey */ FieldIntegratorNull<T>,
  /* Ez */ FieldIntegratorNull<T>,
  /* Hx */ FieldIntegratorNull<T>,
  /* Hy */ FieldIntegratorNull<T>,
  /* Hz */ FieldIntegratorNull<T>
>;

template<typename T>
using EM1D = TypeList<
  /* Ex */ FieldIntegratorNull<T>,
  /* Ey */ FieldIntegratorNull<T>,
  /* Ez */ FieldIntegrator1D<T, FieldUpdate<Derivative::DX, Derivative::NoOp, false, size_t>>,
  /* Hx */ FieldIntegratorNull<T>,
  /* Hy */ FieldIntegrator1D<T, FieldUpdate<Derivative::DX, Derivative::NoOp, true, size_t>>,
  /* Hz */ FieldIntegratorNull<T>
>;

template<typename T>
using EMTM = TypeList<
  /* Ex */ FieldIntegratorNull<T>,
  /* Ey */ FieldIntegratorNull<T>,
  /* Ez */ FieldIntegrator2D<T, FieldUpdate<Derivative::DX, Derivative::DY, false, size_t, size_t>>,
  /* Hx */ FieldIntegrator2D<T, FieldUpdate<Derivative::NoOp, Derivative::DY, true, size_t, size_t>>,
  /* Hy */ FieldIntegrator2D<T, FieldUpdate<Derivative::DX, Derivative::NoOp, true, size_t, size_t>>,
  /* Hz */ FieldIntegratorNull<T>
>;

template<typename T>
using EMTE = TypeList<
  /* Ex */ FieldIntegrator2D<T, FieldUpdate<Derivative::DY, Derivative::NoOp, false, size_t, size_t>>,
  /* Ey */ FieldIntegrator2D<T, FieldUpdate<Derivative::NoOp, Derivative::DX, false, size_t, size_t>>,
  /* Ez */ FieldIntegratorNull<T>,
  /* Hx */ FieldIntegratorNull<T>,
  /* Hy */ FieldIntegratorNull<T>,
  /* Hz */ FieldIntegrator2D<T, FieldUpdate<Derivative::DY, Derivative::DX, true, size_t, size_t>>
>;

template<typename T>
using EM3D = TypeList<
  /* Ex */ FieldIntegrator3D<T, FieldUpdate<Derivative::DY, Derivative::DZ, false, size_t, size_t, size_t>>,
  /* Ey */ FieldIntegrator3D<T, FieldUpdate<Derivative::DZ, Derivative::DX, false, size_t, size_t, size_t>>,
  /* Ez */ FieldIntegrator3D<T, FieldUpdate<Derivative::DX, Derivative::DY, false, size_t, size_t, size_t>>,
  /* Hx */ FieldIntegrator3D<T, FieldUpdate<Derivative::DY, Derivative::DZ, true, size_t, size_t, size_t>>,
  /* Hy */ FieldIntegrator3D<T, FieldUpdate<Derivative::DZ, Derivative::DX, true, size_t, size_t, size_t>>,
  /* Hz */ FieldIntegrator3D<T, FieldUpdate<Derivative::DX, Derivative::DY, true, size_t, size_t, size_t>>
>;

//=================== BCData Definitions ========================
//===============================================================
template<typename T>
using bcdataNone = BCData<
  EmptyArray1D<T>,
  disabled<EmptyArray1D<T>>, // Ex
  disabled<EmptyArray1D<T>>, // Ey
  disabled<EmptyArray1D<T>>, // Ez
  disabled<EmptyArray1D<T>>, // Hx
  disabled<EmptyArray1D<T>>, // Hy
  disabled<EmptyArray1D<T>>  // Hz
>;

template<typename T>
using bcdata1d = BCData<
  Array1D<T>,
  disabled<EmptyArray1D<T>>, // Ex
  disabled<EmptyArray1D<T>>, // Ey
  enabled<Array1D<T>>,       // Ez
  disabled<EmptyArray1D<T>>, // Hx
  enabled<Array1D<T>>,       // Hy
  disabled<EmptyArray1D<T>>  // Hz
>;

template<typename T>
using bcdataTM = BCData<
  Array2D<T>,
  disabled<EmptyArray2D<T>>, // Ex
  disabled<EmptyArray2D<T>>, // Ey
  enabled<Array2D<T>>,       // Ez
  enabled<Array2D<T>>,       // Hx
  enabled<Array2D<T>>,       // Hy
  disabled<EmptyArray2D<T>>  // Hz
>;

template<typename T>
using bcdataTE = BCData<
  Array2D<T>,
  enabled<Array2D<T>>,       // Ex
  enabled<Array2D<T>>,       // Ey
  disabled<EmptyArray2D<T>>, // Ez
  disabled<EmptyArray2D<T>>, // Hx
  disabled<EmptyArray2D<T>>, // Hy
  enabled<Array2D<T>>        // Hz
>;

template<typename T>
using bcdata3D = BCData<
  Array3D<T>,
  enabled<Array3D<T>>, // Ex
  enabled<Array3D<T>>, // Ey
  enabled<Array3D<T>>, // Ez
  enabled<Array3D<T>>, // Hx
  enabled<Array3D<T>>, // Hy
  enabled<Array3D<T>>  // Hz
>;

//=================== Boundary Conditions Definitions ========================
//============================================================================


// enum class EMFace { X, Y, Z};
// enum class EMSide { Lo, Hi};

template<EMFace Face, EMSide Side>
struct Boundary {
  static constexpr EMFace face = Face;
  static constexpr EMSide side = Side;

  template<typename T>
  static IntegratorOffsets get_offsets(const T&, size_t);
};

template <>
template <typename T>
IntegratorOffsets Boundary<EMFace::X, EMSide::Lo>::get_offsets(const T& f, const size_t depth) {
  return {0, depth, 0, f.ny(), 0, f.nz()};
}

template <>
template <typename T>
IntegratorOffsets Boundary<EMFace::X, EMSide::Hi>::get_offsets(const T& f, const size_t depth) {
  return {f.nx() - depth, f.nx(), 0, f.ny(), 0, f.nz()};
}

template <>
template <typename T>
IntegratorOffsets Boundary<EMFace::Y, EMSide::Lo>::get_offsets(const T& f, const size_t depth) {
  return {0, f.nx(), 0, depth, 0, f.nz()};
}

template <>
template <typename T>
IntegratorOffsets Boundary<EMFace::Y, EMSide::Hi>::get_offsets(const T& f, const size_t depth) {
  return {0, f.nx(), f.ny() - depth, f.ny(), 0, f.nz()};
}

template <>
template <typename T>
IntegratorOffsets Boundary<EMFace::Z, EMSide::Lo>::get_offsets(const T& f, const size_t depth) {
  return {0, f.nx(), 0, f.ny(), 0, depth};
}

template <>
template <typename T>
IntegratorOffsets Boundary<EMFace::Z, EMSide::Hi>::get_offsets(const T& f, const size_t depth) {
  return {0, f.nx(), 0, f.ny(), f.nz() - depth, f.nz()};
}

using XLo = Boundary<EMFace::X, EMSide::Lo>;
using YLo = Boundary<EMFace::Y, EMSide::Lo>;
using ZLo = Boundary<EMFace::Z, EMSide::Lo>;
using XHi = Boundary<EMFace::X, EMSide::Hi>;
using YHi = Boundary<EMFace::Y, EMSide::Hi>;
using ZHi = Boundary<EMFace::Z, EMSide::Hi>;

template<typename T, EMFace Face>
using Periodic1D = PeriodicBC<T, Face, nHalo, size_t>;

template<typename T, EMFace Face>
using Periodic2D = PeriodicBC<T, Face, nHalo, size_t, size_t>;

template<typename T, EMFace Face>
using Periodic3D = PeriodicBC<T, Face, nHalo, size_t, size_t, size_t>;

template<typename Boundary, typename Ex, typename Ey, typename Ez, typename Hx, typename Hy, typename Hz>
struct BCApplicator {
  static constexpr EMFace face = Boundary::face;

  static void applyE(auto& emdata, auto& bcdata) {
    static const auto ex_offsets = Boundary::get_offsets(emdata.Ex, Ex::bc_depth);
    static const auto ey_offsets = Boundary::get_offsets(emdata.Ey, Ey::bc_depth);
    static const auto ez_offsets = Boundary::get_offsets(emdata.Ez, Ez::bc_depth);
    // DBG(ez_offsets.x0, ez_offsets.x1, ez_offsets.y0, ez_offsets.y1, ez_offsets.z0, ez_offsets.z1);

    Ex::apply(emdata.Ex, emdata.Hz, emdata.Hy, emdata.Cexh, bcdata.psiEx, bcdata.bEx, bcdata.cEx, ex_offsets);
    Ey::apply(emdata.Ey, emdata.Hx, emdata.Hz, emdata.Ceyh, bcdata.psiEy, bcdata.bEy, bcdata.cEy, ey_offsets);
    Ez::apply(emdata.Ez, emdata.Hy, emdata.Hx, emdata.Cezh, bcdata.psiEz, bcdata.bEz, bcdata.cEz, ez_offsets);
  }

  static void applyH(auto& emdata, auto& bcdata) {
    static const auto hx_offsets = Boundary::get_offsets(emdata.Hx, Hx::bc_depth);
    static const auto hy_offsets = Boundary::get_offsets(emdata.Hy, Hy::bc_depth);
    static const auto hz_offsets = Boundary::get_offsets(emdata.Hz, Hz::bc_depth);
    // DBG(hy_offsets.x0, hy_offsets.x1, hy_offsets.y0, hy_offsets.y1, hy_offsets.z0, hy_offsets.z1);

    Hx::apply(emdata.Hx, emdata.Ey, emdata.Ez, emdata.Chxe, bcdata.psiHx, bcdata.bHx, bcdata.cHx, hx_offsets);
    Hy::apply(emdata.Hy, emdata.Ez, emdata.Ex, emdata.Chye, bcdata.psiHy, bcdata.bHy, bcdata.cHy, hy_offsets);
    Hz::apply(emdata.Hz, emdata.Ex, emdata.Ey, emdata.Chze, bcdata.psiHz, bcdata.bHz, bcdata.cHz, hz_offsets);
  }
};

// template<typename BC>
// requires std::same_as<BC, YLo>
// struct BCApplicator {};
//
// template<typename BC>
// requires std::same_as<BC, ZLo>
// struct BCApplicator {};

#endif //EM_DEFINITIONS_H
