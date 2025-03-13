#ifndef EM_BOUNDARIES_HPP
#define EM_BOUNDARIES_HPP

#include "em_params.hpp"
#include "array.hpp"
#include "bc_functors.hpp"

#include <array>
#include <vector>

namespace tf::electromagnetics {
  //  struct PeriodicData {
  //    using offset_t = std::array<std::size_t, 6>;
  //
  //    explicit PeriodicData(std::size_t nInterior_, std::size_t hi_idx_, offset_t offsets_)
  //    : nInterior(nInterior_), hi_idx(hi_idx_), offsets(offsets_)
  //    {}
  //
  //    std::size_t nInterior;
  //    std::size_t hi_idx;
  //    offset_t offsets;
  //  };

  struct PMLData {
    using offset_t = std::array<std::size_t, 6>;
    using coeffs_t = std::array<double, PMLDepth>;

    explicit PMLData() = delete;
    PMLData(std::size_t, std::size_t, std::size_t, const offset_t&, const offset_t&, bool=false);

    void init_coefficients(bool);
    static std::vector<double> calculate_sigma(const std::vector<double>&, double);
    static std::vector<double> calculate_alpha(const std::vector<double>&);
    static void calculate_coeffs(coeffs_t&, coeffs_t&, const std::vector<double>&, const std::vector<double>&);

    Array3D<double> psiE;
    Array3D<double> psiH;
    offset_t offsetsE;
    offset_t offsetsH;
    coeffs_t bE{};
    coeffs_t cE{};
    coeffs_t bH{};
    coeffs_t cH{};
  };


  struct BCData {
    using Ey_x0_bc = BCIntegrator<double, PMLFunctor<forward_dx, false, true, true>>; // todo: is forward correct?
    using Ey_x1_bc = BCIntegrator<double, PMLFunctor<forward_dx, true, true, true>>;

    using Hy_x0_bc = BCIntegrator<double, PMLFunctor<backward_dx, false, false, false>>; // todo: is forward correct?
    using Hy_x1_bc = BCIntegrator<double, PMLFunctor<backward_dx, true, false, false>>;

    BCData(std::size_t, std::size_t, std::size_t);

    PMLData x0;
    PMLData y0;
    PMLData z0;
    PMLData x1;
    PMLData y1;
    PMLData z1;
  };
}

#endif //EM_BOUNDARIES_HPP
