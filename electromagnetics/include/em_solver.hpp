#ifndef EM_SOLVER_HPP
#define EM_SOLVER_HPP

#include "em_params.hpp"
#include "em_data.hpp"
#include "update_functors.hpp"
#include "bc_functors.hpp"
#include "diff_operators.hpp"

namespace tf::electromagnetics {
  struct EMSolver {
    using ex_func = FieldIntegrator<ex_t, ExplicitUpdateFunctor<backward_dy, backward_dz>>;
    using ey_func = FieldIntegrator<ey_t, ExplicitUpdateFunctor<backward_dz, backward_dx>>;
    using ez_func = FieldIntegrator<ez_t, ExplicitUpdateFunctor<backward_dx, backward_dy>>;
    using hx_func = FieldIntegrator<hx_t, ExplicitUpdateFunctor<forward_dz, forward_dy>>;
    using hy_func = FieldIntegrator<hy_t, ExplicitUpdateFunctor<forward_dx, forward_dz>>;
    using hz_func = FieldIntegrator<hz_t, ExplicitUpdateFunctor<forward_dy, forward_dx>>;

    explicit EMSolver() = delete;
    explicit EMSolver(std::size_t, std::size_t, std::size_t, double, double);

    void advance(double);
    void updateE();
    void updateH();
    void apply_srcs(double) const;

    EMData emdata;
    ex_func ex_update{};
    ey_func ey_update{};
    ez_func ez_update{};
    hx_func hx_update{};
    hy_func hy_update{};
    hz_func hz_update{};
  };
} // end namespace tf::electromagnetics

#endif //EM_SOLVER_HPP
