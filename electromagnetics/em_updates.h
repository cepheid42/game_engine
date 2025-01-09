//
// Created by cepheid on 12/17/24.
//

#ifndef EM_UPDATES_H
#define EM_UPDATES_H

#include "electromagnetics.param"
#include "aydenstuff/array.h"
#include "em_traits.h"
#include "curl_operators.h"
#include "offsets.h"

namespace tf::electromagnetics
{
  using namespace tf::electromagnetics::traits;

  //=================== Field Functors ========================
  //===========================================================
  template<Derivative D1, Derivative D2, bool Forward, typename... IDXS>
  struct FieldUpdate {
    using CurlA = Diff<D1, Forward>;
    using CurlB = Diff<D2, Forward>;

    static void apply(auto& f, const auto& d1, const auto& d2, const auto& j, const auto& c_f, const auto& c_d1, const auto& c_d2, const auto& c_j, IDXS... idxs) {
      const auto    self = c_f(idxs...) * f(idxs...);
      const auto   diff1 = c_d1(idxs...) * CurlA::apply(d1, idxs...);
      const auto   diff2 = c_d2(idxs...) * CurlB::apply(d2, idxs...);
      const auto    diff = diff1 - diff2;
      const auto current = c_j(idxs...) * j(idxs...);
      f(idxs...) = self + diff - current;
    }
  };

  template<typename T, typename UpdateFunctor>
  struct FieldIntegrator3D {
    using value_t = typename T::value_t;
    using dimension_t = typename T::dimension_t;
    using array_t = tf::types::Array3D<value_t>;
    using update_func = UpdateFunctor;
    using offset_t = tf::electromagnetics::types::IntegratorOffsets;

    static auto apply(auto& f, const auto& d1, const auto& d2, const auto& js, const auto& c_f, const auto& c_d1, const auto& c_d2, const auto& c_src, const offset_t& o) {
#pragma omp parallel for collapse(3) num_threads(NTHREADS)
      for (size_t i = o.x0; i < f.nx() - o.x1; ++i) {
        for (size_t j = o.y0; j < f.ny() - o.y1; ++j) {
          for (size_t k = o.z0; k < f.nz() - o.z1; ++k) {
            update_func::apply(f, d1, d2, js, c_f, c_d1, c_d2, c_src, i, j, k);
          }
        }
      }
    }
  };

  template<typename Array>
  struct FieldIntegratorNull {
    using value_t = typename Array::value_t;
    using dimension_t = typename Array::dimension_t;
    using array_t = tf::types::EmptyArray<Array>;

    using offset_t = tf::electromagnetics::types::IntegratorOffsets;
    // static constexpr void apply(const auto&...) {}
    static constexpr void apply(const auto&, const auto&, const auto&, const auto&, const auto&, const auto&, const auto&, const auto&, const offset_t&) {}
  };
} // end namespace tf::electromagnetics
#endif //EM_UPDATES_H
