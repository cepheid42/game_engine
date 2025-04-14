#ifndef EM_SOLVER_HPP
#define EM_SOLVER_HPP

#include "program_params.hpp"
#include "em_data.hpp"
#include "bc_data.hpp"
#include "update_functors.hpp"
#include "bc_functors.hpp"
#include "diff_operators.hpp"

namespace tf::electromagnetics {
  struct EMSolver {
    using ex_func = FieldIntegrator<compute_t, ExplicitUpdateFunctor<backward_dy, backward_dz>>;
    using ey_func = FieldIntegrator<compute_t, ExplicitUpdateFunctor<backward_dz, backward_dx>>;
    using ez_func = FieldIntegrator<compute_t, ExplicitUpdateFunctor<backward_dx, backward_dy>>;
    using hx_func = FieldIntegrator<compute_t, ExplicitUpdateFunctor<forward_dz, forward_dy>>;
    using hy_func = FieldIntegrator<compute_t, ExplicitUpdateFunctor<forward_dx, forward_dz>>;
    using hz_func = FieldIntegrator<compute_t, ExplicitUpdateFunctor<forward_dy, forward_dx>>;

    // // X-Faces
    // using Ey_x0_bc = BCIntegrator<void>;
    // using Ey_x1_bc = BCIntegrator<void>;
    // using Ez_x0_bc = BCIntegrator<void>;
    // using Ez_x1_bc = BCIntegrator<void>;
    //
    // using Hy_x0_bc = BCIntegrator<PeriodicFunctor<EMFace::X>>;
    // using Hy_x1_bc = BCIntegrator<void>;
    // using Hz_x0_bc = BCIntegrator<PeriodicFunctor<EMFace::X>>;
    // using Hz_x1_bc = BCIntegrator<void>;
    //
    // // Y-Faces
    // using Ex_y0_bc = BCIntegrator<void>;
    // using Ex_y1_bc = BCIntegrator<void>;
    // using Ez_y0_bc = BCIntegrator<void>;
    // using Ez_y1_bc = BCIntegrator<void>;
    //
    // using Hx_y0_bc = BCIntegrator<PeriodicFunctor<EMFace::Y>>;
    // using Hx_y1_bc = BCIntegrator<void>;
    // using Hz_y0_bc = BCIntegrator<PeriodicFunctor<EMFace::Y>>;
    // using Hz_y1_bc = BCIntegrator<void>;
    //
    // // Z-Faces
    // using Ex_z0_bc = BCIntegrator<void>;
    // using Ex_z1_bc = BCIntegrator<void>;
    // using Ey_z0_bc = BCIntegrator<void>;
    // using Ey_z1_bc = BCIntegrator<void>;
    //
    // using Hx_z0_bc = BCIntegrator<PeriodicFunctor<EMFace::Z>>;
    // using Hx_z1_bc = BCIntegrator<void>;
    // using Hy_z0_bc = BCIntegrator<PeriodicFunctor<EMFace::Z>>;
    // using Hy_z1_bc = BCIntegrator<void>;
    
    // X-Faces
    using Ex_x0_bc = BCIntegrator<void>;
    using Ex_x1_bc = BCIntegrator<void>;
    using Ey_x0_bc = BCIntegrator<PMLFunctor<backward_dx, false,  true>>;
    using Ey_x1_bc = BCIntegrator<PMLFunctor<backward_dx,  true,  true>>;
    using Ez_x0_bc = BCIntegrator<PMLFunctor<backward_dx, false, false>>;
    using Ez_x1_bc = BCIntegrator<PMLFunctor<backward_dx,  true, false>>;

    using Hx_x0_bc = BCIntegrator<void>;
    using Hx_x1_bc = BCIntegrator<void>;
    using Hy_x0_bc = BCIntegrator<PMLFunctor<forward_dx, false, false>>;
    using Hy_x1_bc = BCIntegrator<PMLFunctor<forward_dx,  true, false>>;
    using Hz_x0_bc = BCIntegrator<PMLFunctor<forward_dx, false,  true>>;
    using Hz_x1_bc = BCIntegrator<PMLFunctor<forward_dx,  true,  true>>;

    // Y-Faces
    using Ex_y0_bc = BCIntegrator<void>;
    using Ex_y1_bc = BCIntegrator<void>;
    using Ey_y0_bc = BCIntegrator<void>;
    using Ey_y1_bc = BCIntegrator<void>;
    using Ez_y0_bc = BCIntegrator<void>;
    using Ez_y1_bc = BCIntegrator<void>;

    using Hx_y0_bc = BCIntegrator<void>;
    using Hx_y1_bc = BCIntegrator<void>;
    using Hy_y0_bc = BCIntegrator<void>;
    using Hy_y1_bc = BCIntegrator<void>;
    using Hz_y0_bc = BCIntegrator<void>;
    using Hz_y1_bc = BCIntegrator<void>;

    // Z-Faces
    using Ex_z0_bc = BCIntegrator<PMLFunctor<backward_dz, false,  true>>;
    using Ex_z1_bc = BCIntegrator<PMLFunctor<backward_dz,  true,  true>>;
    using Ey_z0_bc = BCIntegrator<PMLFunctor<backward_dz, false, false>>;
    using Ey_z1_bc = BCIntegrator<PMLFunctor<backward_dz,  true, false>>;
    using Ez_z0_bc = BCIntegrator<void>;
    using Ez_z1_bc = BCIntegrator<void>;

    using Hx_z0_bc = BCIntegrator<PMLFunctor<forward_dz, false, false>>;
    using Hx_z1_bc = BCIntegrator<PMLFunctor<forward_dz,  true, false>>;
    using Hy_z0_bc = BCIntegrator<PMLFunctor<forward_dz, false,  true>>;
    using Hy_z1_bc = BCIntegrator<PMLFunctor<forward_dz,  true,  true>>;
    using Hz_z0_bc = BCIntegrator<void>;
    using Hz_z1_bc = BCIntegrator<void>;


    explicit EMSolver() = delete;
    explicit EMSolver(std::size_t, std::size_t, std::size_t, compute_t, compute_t);

    void advance(compute_t);
    void updateE();
    void updateH();
    void updateEBCs();
    void updateHBCs();
    void apply_srcs(compute_t) const;

    EMData emdata;
    BCData<BoundaryType> bcdata;
    
    ex_func ex_update{};
    ey_func ey_update{};
    ez_func ez_update{};
    hx_func hx_update{};
    hy_func hy_update{};
    hz_func hz_update{};

    Ey_x0_bc Ey_x0{};
    Ey_x1_bc Ey_x1{};
    Ez_x0_bc Ez_x0{};
    Ez_x1_bc Ez_x1{};

    Hy_x0_bc Hy_x0{};
    Hy_x1_bc Hy_x1{};
    Hz_x0_bc Hz_x0{};
    Hz_x1_bc Hz_x1{};
    
    Ex_y0_bc Ex_y0{};
    Ex_y1_bc Ex_y1{};
    Ez_y0_bc Ez_y0{};
    Ez_y1_bc Ez_y1{};

    Hx_y0_bc Hx_y0{};
    Hx_y1_bc Hx_y1{};
    Hz_y0_bc Hz_y0{};
    Hz_y1_bc Hz_y1{};
    
    Ex_z0_bc Ex_z0{};
    Ex_z1_bc Ex_z1{};
    Ey_z0_bc Ey_z0{};
    Ey_z1_bc Ey_z1{};

    Hx_z0_bc Hx_z0{};
    Hx_z1_bc Hx_z1{};
    Hy_z0_bc Hy_z0{};
    Hy_z1_bc Hy_z1{};
  };
} // end namespace tf::electromagnetics

#endif //EM_SOLVER_HPP
