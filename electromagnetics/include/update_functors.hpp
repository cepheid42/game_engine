#ifndef EM_UPDATES_HPP
#define EM_UPDATES_HPP

#include <array>
#include <concepts>

namespace tf::electromagnetics {
  template<typename T, typename UpdateFunc>
  struct FieldIntegrator {
    using offset_t = std::array<std::size_t, 6>;

    // todo: make sure variadic args here don't slow anything down like that one time
    static constexpr void operator()(auto&...)
    requires std::same_as<T, void>
    {}

    static void operator()(auto& f, const auto& d1, const auto& d2, const auto& src,
                           const auto& c_f, const auto& c_d1, const auto& c_d2, const auto& c_src,
                           const offset_t& offsets)
    {
      const auto& [x0, x1, y0, y1, z0, z1] = offsets;
#pragma omp parallel for simd collapse(3) num_threads(nThreads)
      for (std::size_t i = x0; i < f.nx() - x1; ++i) {
        for (std::size_t j = y0; j < f.ny() - y1; ++j) {
          for (std::size_t k = z0; k < f.nz() - z1; ++k) {
            UpdateFunc::apply(f, d1, d2, src, c_f, c_d1, c_d2, c_src, i, j, k);
          } // end for k
        } // end for j
      } // end for i
    } // end operator()
  }; // end struct FieldIntegrator

  template<typename CurlA, typename CurlB>
  struct ExplicitUpdateFunctor {
    static void apply(auto& f, const auto& d1, const auto& d2, const auto& src,
                      const auto& c_f, const auto& c_d1, const auto& c_d2, const auto& c_src,
                      const std::size_t i, const std::size_t j, const std::size_t k)
    {
      const auto    self = c_f(i, j, k) * f(i, j, k);
      const auto   diff1 = c_d1(i, j, k) * CurlA::apply(d1, i, j, k);
      const auto   diff2 = c_d2(i, j, k) * CurlB::apply(d2, i, j, k);
      const auto current = c_src(i, j, k) * src(i, j, k);
      f(i, j, k) = self + (diff1 - diff2) - current;
    } // end apply
  }; // end struct ExplicitUpdateFunctor

//  struct ImplicitUpdateFunctor {
//    static void apply() {
//      // todo : tridiagonal solve here
//    }
//  };
//
//  template<typename CurlA, typename CurlB>
//  struct ExplicitFieldUpdate {
//    static void apply(auto& f, const auto& d1, const auto& d2, const auto& src,
//                      const auto& c_f, const auto& c_d1, const auto& c_d2, const auto& c_src,
//                      const std::size_t i, const std::size_t j, const std::size_t k)
//    {
//      ExplicitUpdateFunctor<CurlA, CurlB>::apply(f, d1, d2, src, c_f, c_d1, c_d2, c_src, i, j, k);
//    }
//  };
//
//  template<typename CurlA, typename CurlB>
//  struct ImplicitFieldUpdate {
//    static void apply(auto& f, auto& f_aux, const auto& d1, const auto& d2, const auto& src,
//                      const auto& c_f, const auto& c_d1, const auto& c_d2, const auto& c_src,
//                      const std::size_t i, const std::size_t j, const std::size_t k)
//    {
//      // RHS
//      ExplicitUpdateFunctor<CurlA, CurlB>::apply(f_aux, d1, d2, src, c_f, c_d1, c_d2, c_src, i, j, k);
//
//      // LHS
//      ImplicitUpdateFunctor::apply();
//
//      // rest of the update
//    }
//  };
}



#endif //EM_UPDATES_HPP
