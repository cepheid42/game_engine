#ifndef EM_UPDATES_HPP
#define EM_UPDATES_HPP

#include <array>

namespace tf::electromagnetics {

template<typename CurlA, typename CurlB>
struct ExplicitUpdateFunctor {
   #pragma omp declare simd notinbranch
   static void apply(auto& f, const auto& d1, const auto& d2, const auto& src,
                     const auto& c_f, const auto& c_d1, const auto& c_d2, const auto& c_src,
                     const std::size_t i, const std::size_t j, const std::size_t k)
   {
      // const auto self    = c_f(i, j, k) * f(i, j, k);
      // const auto diff1   = c_d1(i, j, k) * CurlA::apply(d1, i, j, k);
      // const auto diff2   = c_d2(i, j, k) * CurlB::apply(d2, i, j, k);
      // const auto current = c_src(i, j, k) * src(i, j, k);
      // f(i, j, k)         = self + (diff1 - diff2) - current;
      f(i, j, k) =  c_f(i, j, k) * f(i, j, k)
                 + (c_d1(i, j, k) * CurlA::apply(d1, i, j, k)
                 -  c_d2(i, j, k) * CurlB::apply(d2, i, j, k))
                 -  c_src(i, j, k) * src(i, j, k);
   } // end apply
}; // end struct ExplicitUpdateFunctor


template<typename UpdateFunc>
struct FieldIntegrator {
   using offset_t = std::array<std::size_t, 6>;

   static void apply(auto& f, const auto& d1, const auto& d2, const auto& src,
                     const auto& c_f, const auto& c_d1, const auto& c_d2, const auto& c_src,
                     const offset_t& offsets)
   {
      #pragma omp parallel for simd collapse(3) num_threads(nThreads)
      for (std::size_t i = offsets[0]; i < f.nx() - offsets[1]; ++i) {
         for (std::size_t j = offsets[2]; j < f.ny() - offsets[3]; ++j) {
            for (std::size_t k = offsets[4]; k < f.nz() - offsets[5]; ++k) {
               UpdateFunc::apply(f, d1, d2, src, c_f, c_d1, c_d2, c_src, i, j, k);
            } // end for k
         } // end for j
      } // end for i
   } // end operator()
}; // end struct FieldIntegrator

} // end namespace tf::electromagnetics


#endif //EM_UPDATES_HPP
