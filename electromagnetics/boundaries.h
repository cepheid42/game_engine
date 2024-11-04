//
// Created by cepheid on 10/17/24.
//

#ifndef BOUNDARIES_H
#define BOUNDARIES_H


// #include "aydenstuff/array.h"
// #include <electromagnetics.h>

// #include "bc_data.h"
#include "offsets.h"


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
  using value_t = typename Array::value_t;
  using dimension_t = typename Array::dimension_t;

  static constexpr size_t bc_depth = DEPTH;

  static void apply1D(auto& f1, const size_t i) {
    const auto numInterior = f1.nx() - (2 * bc_depth);
    const auto hi_idx = (f1.nx() - 1) - bc_depth;
    const auto pm = i % numInterior;

    f1[bc_depth - 1 - i] = f1[hi_idx - pm];
    f1[hi_idx + 1 + i] = f1[bc_depth + pm];
  }

  static void apply2D(auto& f1, const size_t i, const size_t j) {
    if constexpr (Face == EMFace::X) {
      const auto numInterior = f1.nx() - (2 * bc_depth);
      const auto hi_idx = (f1.nx() - 1) - bc_depth;
      const auto pm = i % numInterior;

      f1(bc_depth - 1 - i, j) = f1(hi_idx - pm, j);
      f1(hi_idx + 1 + i, j) = f1(bc_depth + pm, j);

    } else if constexpr (Face == EMFace::Y) {
      const auto numInterior = f1.ny() - (2 * bc_depth);
      const auto hi_idx = (f1.ny() - 1) - bc_depth;
      const auto pm = j % numInterior;

      f1(i, bc_depth - 1 - j) = f1(i, hi_idx - pm);
      f1(i, hi_idx + 1 + j) = f1(i, bc_depth + pm);
    }
  }


  static void apply3D(auto& f1, const size_t i, const size_t j, const size_t k) {
    if constexpr (Face == EMFace::X) {
      // const auto numInterior = f1.nx() - (2 * bc_depth);
      const auto hi_idx = (f1.nx() - 1) - bc_depth;
      const auto pm = i % (f1.nx() - (2 * bc_depth)); // i % numInterior

      f1(bc_depth - 1 - i, j, k) = f1(hi_idx - pm, j, k);
      f1(hi_idx + 1 + i, j, k) = f1(bc_depth + pm, j, k);
    } else if constexpr (Face == EMFace::Y) {
      // const auto numInterior = f1.ny() - (2 * bc_depth);
      const auto hi_idx = (f1.ny() - 1) - bc_depth;
      const auto pm = j % (f1.ny() - (2 * bc_depth)); // j % numInterior

      f1(i, bc_depth - 1 - j, k) = f1(i, hi_idx - pm, k);
      f1(i, hi_idx + 1 + j, k) = f1(i, bc_depth + pm, k);
    } else if constexpr (Face == EMFace::Z) {
      // const auto numInterior = f1.nz() - (2 * bc_depth);
      const auto hi_idx = (f1.nz() - 1) - bc_depth;
      const auto pm = k % (f1.nz() - (2 * bc_depth)); // k % numInterior

      f1(i, j, bc_depth - 1 - k) = f1(i, j, hi_idx - pm);
      f1(i, j, hi_idx + 1 + k) = f1(i, j, bc_depth + pm);
    }
  }

  template<typename ARR>
  static void apply(ARR& f1, const auto&, const auto&, const auto&, const auto&, const auto&, const auto&, IDXS... idxs) {
    // DBG("PeriodicBC::apply()");
    if constexpr (dimension_t::value == 1) {
      apply1D(f1, idxs...);
    } else if constexpr (dimension_t::value == 2) {
      apply2D(f1, idxs...);
    } else {
      apply3D(f1, idxs...);
    }
  }

  template<typename ARR>
  requires is_empty_field<ARR, EmptyArray<value_t, dimension_t::value>>
  static void apply(ARR&, const auto&, const auto&, const auto&, const auto&, const auto&, const auto&, IDXS...)
  { DBG("PeriodicBC::apply()::empty"); }
};

template<typename Array, EMSide Side, typename PMLData, size_t DEPTH, typename... IDXS>
struct PmlBC {
  using value_t = typename Array::value_t;
  using dimension_t = typename Array::dimension_t;

  using CurlA = typename PMLData::curl1;
  using CurlB = typename PMLData::curl2;

  static constexpr size_t bc_depth = DEPTH;

  static void apply1D(auto& f1, const auto& d1, const auto&, const auto& c_d, auto& psi, const auto& b, const auto& c, const size_t i) {
    static size_t start;
    if constexpr (Side == EMSide::Hi) {
      start = (f1.nx - 1) - bc_depth;  // need -1?
    } else {
      start = 0;
    }
    const auto ipml = i - start;

    const auto self = b[ipml] * psi[ipml];
    const auto diff1 = CurlA::apply(d1, i);
    const auto diff2 = CurlB::apply(d1, i);
    const auto diff = c[ipml] * (diff1 - diff2);

    psi[ipml] = self + diff;
    f1[i] += c_d[i] * psi[ipml];
  }


  static void apply2D(auto& f1, const auto& d1, const auto& d2, const auto& c_d, auto& psi, const auto& b, const auto& c, const size_t i, const size_t j) {}
  static void apply3D(auto& f1, const auto& d1, const auto& d2, const auto& c_d, auto& psi, const auto& b, const auto& c, const size_t i, const size_t j, const size_t k) {}

  static void apply(auto& f1, const auto& d1, const auto& d2, const auto& c_d, auto& psi, const auto& b, const auto& c, IDXS... idxs) {
    if constexpr (dimension_t::value == 1) {
      apply1D(f1, d1, d2, c_d, psi, b, c, idxs...);
    } else if constexpr (dimension_t::value == 2) {
      apply2D(f1, d1, d2, c_d, psi, b, c, idxs...);
    } else {
      apply3D(f1, d1, d2, c_d, psi, b, c, idxs...);
    }
  }

  template<typename ARR>
  requires is_empty_field<ARR, EmptyArray<value_t, dimension_t::value>>
  static void apply(ARR&, const auto&, const auto&, const auto&, const auto&, const auto&)
  { DBG("PML1D::apply()::empty"); }
};


template<typename T, typename UpdateFunctor>
struct BCIntegrator1D {
  using value_t = typename T::value_t;
  using dimension_t = typename T::dimension_t;
  using array_t = typename T::array_t;
  using update_func = UpdateFunctor;
  using offset_t = IntegratorOffsets;
  static constexpr size_t bc_depth = UpdateFunctor::bc_depth;

  static void apply(auto& f1, const auto& d1, const auto& d2, const auto& c_d, auto& psi, const auto& b, const auto& c, const offset_t& o) {
    if constexpr (!std::same_as<UpdateFunctor, ReflectingBC<array_t>>) {
      DBG("BCIntegrator1D::apply()");
      for (size_t i = o.x0; i < o.x1; ++i) {
        update_func::apply(f1, d1, d2, c_d, psi, b, c, i);
      }
    } else {
      DBG("BCIntegrator1D::apply()::empty");
    }
  }
};

template<typename T, typename UpdateFunctor>
struct BCIntegrator2D {
  using value_t = typename T::value_t;
  using dimension_t = typename T::dimension_t;
  using array_t = Array2D<value_t>;
  using update_func = UpdateFunctor;
  using offset_t = IntegratorOffsets;
  static constexpr size_t bc_depth = UpdateFunctor::bc_depth;

  static void apply(auto& f1, const auto& d1, const auto& d2, const auto& c_d, auto& psi, const auto& b, const auto& c, const offset_t& o) {
    if constexpr (!std::same_as<UpdateFunctor, ReflectingBC<array_t>>) {
      DBG("BCIntegrator2D::apply()");
      // DBG(o.x0, o.x1, o.y0, o.y1);
      for (size_t i = o.x0; i < o.x1; ++i) {
        for (size_t j = o.y0; j < o.y1; ++j) {
          update_func::apply(f1, d1, d2, c_d, psi, b, c, i, j);
        }
      }
    } else {
      DBG("BCIntegrator2D::apply()::empty");
    }
  }
};

template<typename T, typename UpdateFunctor>
struct BCIntegrator3D {
  using value_t = typename T::value_t;
  using dimension_t = typename T::dimension_t;
  using array_t = Array3D<value_t>;
  using update_func = UpdateFunctor;
  using offset_t = IntegratorOffsets;
  static constexpr size_t bc_depth = UpdateFunctor::bc_depth;

  static void apply(auto& f1, const auto& d1, const auto& d2, const auto& c_d, auto& psi, const auto& b, const auto& c, const offset_t& o) {
    if constexpr (!std::same_as<UpdateFunctor, ReflectingBC<array_t>>) {
      DBG("BCIntegrator3D::apply()");
      for (size_t i = o.x0; i < o.x1; ++i) {
        for (size_t j = o.y0; j < o.y1; ++j) {
          for (size_t k = o.z0; k < o.z1; ++k) {
            update_func::apply(f1, d1, d2, c_d, psi, b, c, i, j, k);
          }
        }
      }
    } else {
      DBG("BCIntegrator3D::apply()::empty");
    }
  }
};


#endif //BOUNDARIES_H
