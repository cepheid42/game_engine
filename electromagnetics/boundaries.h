//
// Created by cepheid on 10/17/24.
//

#ifndef BOUNDARIES_H
#define BOUNDARIES_H


#include "offsets.h"

constexpr size_t dPML = 10u;
constexpr size_t nHalo = 2u;


//====== 1D Boundaries ========
//=============================
template<typename Array>
struct ReflectingBC {
  using array_t = Array;
  using value_t = typename Array::value_t;
  using dimension_t = typename Array::dimension_t;
  static constexpr size_t bc_depth = 0;
};

template<typename Array, EMFace Face, size_t DEPTH, typename... IDXS>
struct PeriodicBC {
  // using value_t = typename Array::value_t;
  // using dimension_t = typename Array::dimension_t;

  static constexpr size_t bc_depth = DEPTH;

  // static void apply1D(auto& f1, const size_t i) {
  //   const auto numInterior = f1.nx() - (2 * bc_depth);
  //   const auto hi_idx = (f1.nx() - 1) - bc_depth;
  //   const auto pm = i % numInterior;
  //
  //   f1[bc_depth - 1 - i] = f1[hi_idx - pm];
  //   f1[hi_idx + 1 + i] = f1[bc_depth + pm];
  // }
  //
  // static void apply2D(auto& f1, const size_t i, const size_t j) {
  //   if constexpr (Face == EMFace::X) {
  //     const auto numInterior = f1.nx() - (2 * bc_depth);
  //     const auto hi_idx = (f1.nx() - 1) - bc_depth;
  //     const auto pm = i % numInterior;
  //
  //     f1(bc_depth - 1 - i, j) = f1(hi_idx - pm, j);
  //     f1(hi_idx + 1 + i, j) = f1(bc_depth + pm, j);
  //
  //   } else if constexpr (Face == EMFace::Y) {
  //     const auto numInterior = f1.ny() - (2 * bc_depth);
  //     const auto hi_idx = (f1.ny() - 1) - bc_depth;
  //     const auto pm = j % numInterior;
  //
  //     f1(i, bc_depth - 1 - j) = f1(i, hi_idx - pm);
  //     f1(i, hi_idx + 1 + j) = f1(i, bc_depth + pm);
  //   }
  // }
  //
  //
  // static void apply3D(auto& f1, const size_t i, const size_t j, const size_t k) {
  //   if constexpr (Face == EMFace::X) {
  //     // const auto numInterior = f1.nx() - (2 * bc_depth);
  //     const auto hi_idx = (f1.nx() - 1) - bc_depth;
  //     const auto pm = i % (f1.nx() - (2 * bc_depth)); // i % numInterior
  //
  //     f1(bc_depth - 1 - i, j, k) = f1(hi_idx - pm, j, k);
  //     f1(hi_idx + 1 + i, j, k) = f1(bc_depth + pm, j, k);
  //   } else if constexpr (Face == EMFace::Y) {
  //     // const auto numInterior = f1.ny() - (2 * bc_depth);
  //     const auto hi_idx = (f1.ny() - 1) - bc_depth;
  //     const auto pm = j % (f1.ny() - (2 * bc_depth)); // j % numInterior
  //
  //     f1(i, bc_depth - 1 - j, k) = f1(i, hi_idx - pm, k);
  //     f1(i, hi_idx + 1 + j, k) = f1(i, bc_depth + pm, k);
  //   } else if constexpr (Face == EMFace::Z) {
  //     // const auto numInterior = f1.nz() - (2 * bc_depth);
  //     const auto hi_idx = (f1.nz() - 1) - bc_depth;
  //     const auto pm = k % (f1.nz() - (2 * bc_depth)); // k % numInterior
  //
  //     f1(i, j, bc_depth - 1 - k) = f1(i, j, hi_idx - pm);
  //     f1(i, j, hi_idx + 1 + k) = f1(i, j, bc_depth + pm);
  //   }
  // }

  template<typename ARR>
  static void apply(ARR& f1, const size_t i) {//, const auto&, const auto&, const auto&, const auto&, const auto&, const auto&, IDXS... idxs) {
    DBG("PeriodicBC::apply()");
    // if constexpr (ARR::dimension_t::value == 1) {
    //   apply1D(f1, idxs...);
    // } else if constexpr (ARR::dimension_t::value == 2) {
    //   apply2D(f1, idxs...);
    // } else {
    //   apply3D(f1, idxs...);
    // }
  }

  template<typename ARR>
  requires is_empty_field<ARR, EmptyArray<typename ARR::value_t, ARR::dimension_t::value>>
  static void apply(ARR&, const size_t) //const auto&, const auto&, const auto&, const auto&, const auto&, const auto&, IDXS...)
  { DBG("PeriodicBC::apply()::empty"); }
};

template<typename Array, typename CurlA, typename CurlB, size_t DEPTH, typename... IDXS>
struct PmlBC {
  using array_t = Array;
  using value_t = typename Array::value_t;
  using dimension_t = typename Array::dimension_t;

  static constexpr size_t bc_depth = DEPTH;

  // static void apply1D(auto& f1, const auto& d1, const auto& c_d, auto& psi, const auto& b, const auto& c, const size_t i) {
  //   // static size_t start;
  //   // if constexpr (Side == EMSide::Hi) {
  //   //   start = (f1.nx - 1) - bc_depth;  // need -1?
  //   // } else {
  //   //   start = 0;
  //   // }
  //   const auto ipml = i - start;
  //
  //   const auto self = b[ipml] * psi[ipml];
  //   const auto diff1 = CurlA::apply(d1, i);
  //   const auto diff2 = CurlB::apply(d1, i);
  //   const auto diff = c[ipml] * (diff1 - diff2);
  //
  //   psi[ipml] = self + diff;
  //   f1[i] += c_d[i] * psi[ipml];
  // }


  // static void apply2D(auto& f1, const auto& d1, const auto& d2, const auto& c_d, auto& psi, const auto& b, const auto& c, const size_t i, const size_t j) {}
  // static void apply3D(auto& f1, const auto& d1, const auto& d2, const auto& c_d, auto& psi, const auto& b, const auto& c, const size_t i, const size_t j, const size_t k) {}

  static void apply(auto& f1, size_t i) {//const auto& d1, const auto& d2, const auto& c_d, auto& psi, const auto& b, const auto& c, IDXS... idxs) {
    DBG("PMLBC::apply()");
  }

  template<typename ARR>
  requires is_empty_field<ARR, EmptyArray<value_t, dimension_t::value>>
  static void apply(ARR&, const size_t) // const auto&, const auto&, const auto&, const auto&, const auto&)
  { DBG("PMLBC::apply()::empty"); }
};


template<typename UpdateFunctor>
struct BCIntegrator1D {
  using update_t = UpdateFunctor;
  using array_t = typename update_t::array_t;
  using value_t = typename update_t::value_t;
  using dimension_t = typename update_t::dimension_t;

  static constexpr size_t bc_depth = update_t::bc_depth;

  static void apply(auto& f1) {//, const auto& d1, const auto& d2, const auto& c_d, auto& psi, const auto& b, const auto& c, const offset_t& o) {
    if constexpr (!std::same_as<UpdateFunctor, ReflectingBC<array_t>>) {
      DBG("BCIntegrator1D::apply()");
      for (size_t i = f1.offsets.x0; i < f1.offsets.x1; ++i) {
        update_t::apply(f1, i);
      }
    } else {
      DBG("BCIntegrator1D::apply()::empty");
    }
  }
};

template<typename UpdateFunctor>
struct BCIntegrator2D {
  using update_t = UpdateFunctor;
  using value_t = typename update_t::value_t;
  using dimension_t = typename update_t::dimension_t;
  using array_t = Array2D<value_t>;

  static constexpr size_t bc_depth = UpdateFunctor::bc_depth;

  static void apply(auto& f1) {//, const auto& d1, const auto& d2, const auto& c_d, auto& psi, const auto& b, const auto& c, const offset_t& o) {
    if constexpr (!std::same_as<UpdateFunctor, ReflectingBC<array_t>>) {
      DBG("BCIntegrator2D::apply()");
      for (size_t i = f1.offsets.x0; i < f1.offsets.x1; ++i) {
        for (size_t j = f1.offsets.y0; j < f1.offsets.y1; ++j) {
          update_t::apply(f1, i, j);
        }
      }
    } else {
      DBG("BCIntegrator2D::apply()::empty");
    }
  }
};

template<typename UpdateFunctor>
struct BCIntegrator3D {
  using update_t = UpdateFunctor;
  using value_t = typename update_t::value_t;
  using dimension_t = typename update_t::dimension_t;
  using array_t = Array3D<value_t>;

  static constexpr size_t bc_depth = UpdateFunctor::bc_depth;

  static void apply(auto& f1) { //, const auto& d1, const auto& d2, const auto& c_d, auto& psi, const auto& b, const auto& c, const offset_t& o) {
    if constexpr (!std::same_as<UpdateFunctor, ReflectingBC<array_t>>) {
      DBG("BCIntegrator3D::apply()");
      for (size_t i = f1.offsets.x0; i < f1.offsets.x1; ++i) {
        for (size_t j = f1.offsets.y0; j < f1.offsets.y1; ++j) {
          for (size_t k = f1.offsets.z0; k < f1.offsets.z1; ++k) {
            update_t::apply(f1, i, j, k);
          }
        }
      }
    } else {
      DBG("BCIntegrator3D::apply()::empty");
    }
  }
};

template<typename Integrator>
struct PMLData {
  using array_t = typename Integrator::array_t;
  using value_t = typename Integrator::value_t;
  using dimension_t = typename Integrator::dimension_t;
  using offset_t = IntegratorOffsets;
  using integrator_t = Integrator;

  static constexpr size_t bc_depth = Integrator::bc_depth;

  explicit PMLData(const offset_t& o)
  : offsets{o}, psi{o.x1 - o.x0}, b{o.y1 - o.y0}, c{o.z1 - o.z0}
  {}

  offset_t offsets;
  array_t psi;
  array_t b;
  array_t c;
};

template<typename Integrator>
struct PeriodicData {
  using array_t = typename Integrator::array_t;
  using value_t = typename Integrator::value_t;
  using dimension_t = typename Integrator::dimension_t;
  using offset_t = IntegratorOffsets;

  explicit PeriodicData(const offset_t& o)
  : offsets{o}
  {}

  IntegratorOffsets offsets;
};

// template<typename Array, bool Forward, typename... IDXS>
// using EzHy_XFace_PML = PmlBC<
//   Array,
//   curl<Derivative::NoOp, Forward, IDXS...>,
//   curl<Derivative::DY, Forward, IDXS...>,
//   dPML,
//   IDXS...
// >;
//
// template<typename Integrator>
// using EzBC_X = PMLData<Integrator>;

#endif //BOUNDARIES_H
