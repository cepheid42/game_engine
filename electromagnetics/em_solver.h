//
// Created by cepheid on 9/23/24.
//

#ifndef EM_SOLVER_H
#define EM_SOLVER_H

#include <cassert>
#include <iostream>

#include "./core/typelist.h"
#include "em_data.h"

template<typename T>
void print_type() {
  std::cout << __PRETTY_FUNCTION__ << std::endl;
}

//============== Functor Classes =============
//============================================
struct ForwardDiff1D {
  static auto apply(const auto& d, const std::size_t i) {
    std::cout << "ForwardDiff1D::apply()" << std::endl;
    return 0.0;
    // return d[i + 1] - d[i];
  }
};

struct BackwardDiff1D {
  static auto apply(const auto& d, const std::size_t i) {
    std::cout << "BackwardDiff1D::apply()" << std::endl;
    return 0.0;
    // return d[i] - d[i - 1];
  }
};

//============== Functor Classes =============
//============================================
template<typename emdata_t, typename DIFF>
struct UpdateFunctor1D {
  using value_t = typename emdata_t::value_t;
  using array_t = typename emdata_t::array_t;

  static value_t apply(array_t& F, const array_t& D1, const array_t& J, const array_t& c_F, const array_t& c_D, const array_t& c_J, const std::size_t i) {
    std::cout << "UpdateFunctor1D::apply()" << std::endl;
    DIFF::apply(D1, i);
    return 0.0;
    // const auto    self = c_F[i] * F[i];
    // const auto    diff = c_D[i] * DIFF::apply(D1);
    // const auto current = c_J[i] * J[i];
    // return self + diff - current;
  }

  static value_t apply(array_t& F, const array_t& D1, const array_t& c_F, const array_t& c_D, const std::size_t i) {
    std::cout << "UpdateFunctor1D::apply()" << std::endl;
    DIFF::apply(D1, i);
    return 0.0;
    // const auto    self = c_F[i] * F[i];
    // const auto    diff = c_D[i] * DIFF::apply(D1);
    // const auto current = c_J[i] * J[i];
    // return self + diff - current;
  }
};

//============== Integrator Classes =============
//===============================================
template<typename emdata_t, typename DIFF>
struct Integrator1D {
  using value_t = typename emdata_t::value_t;
  using array_t = typename emdata_t::array_t;

  static void advance(array_t& F, const array_t& D1, const array_t& J, const array_t& c_F, const array_t& c_D, const array_t& c_J) {
    for (std::size_t i = 0; i < 2; i++) {
      UpdateFunctor1D<emdata_t, DIFF>::apply(F, D1, J, c_F, c_D, c_J, i);
    }
  }

  static void advance(array_t& F, const array_t& D1, const array_t& c_F, const array_t& c_D) {
    for (std::size_t i = 0; i < 2; i++) {
      UpdateFunctor1D<emdata_t, DIFF>::apply(F, D1, c_F, c_D, i);
    }
  }
};

struct Integrator2D {
  static void advance() { std::cout << "Integrator2D::advance()" << std::endl; }
};

struct Integrator3D {
  static void advance() { std::cout << "Integrator3D::advance()" << std::endl; }
};

struct NullIntegrator {
  static void advance() { std::cout << "NullIntegrator::advance()" << std::endl; }
};


struct DefaultSolver {
  using EIX = NullIntegrator;
  using EIY = NullIntegrator;
  using EIZ = NullIntegrator;
  using HIX = NullIntegrator;
  using HIY = NullIntegrator;
  using HIZ = NullIntegrator;
};

template<typename Base, int D>
struct Discriminator : Base {};

template<typename ExSetter, typename EySetter, typename EzSetter, typename HxSetter, typename HySetter, typename HzSetter>
struct SolverSelector
: Discriminator<ExSetter, 1>, Discriminator<EySetter, 2>, Discriminator<EzSetter, 3>,
  Discriminator<HxSetter, 4>, Discriminator<HySetter, 5>, Discriminator<HzSetter, 6>
{};

struct DefaultSolverArgs : virtual DefaultSolver {};


template<typename Integrator>
struct ExIntegrator_is : virtual DefaultSolver {
  using EIX = Integrator;
};

template<typename Integrator>
struct EyIntegrator_is : virtual DefaultSolver {
  using EIY = Integrator;
};

template<typename Integrator>
struct EzIntegrator_is : virtual DefaultSolver {
  using EIZ = Integrator;
};

template<typename Integrator>
struct HxIntegrator_is : virtual DefaultSolver {
  using HIX = Integrator;
};

template<typename Integrator>
struct HyIntegrator_is : virtual DefaultSolver {
  using HIY = Integrator;
};

template<typename Integrator>
struct HzIntegrator_is : virtual DefaultSolver {
  using HIZ = Integrator;
};

//============== Electromagnetics Superclass =============
//========================================================
template<typename Solver>
struct EM {
  using emdata_t = typename Solver::emdata_t;

  using EIX = typename Solver::Integrators::EIX;
  using EIY = typename Solver::Integrators::EIY;
  using EIZ = typename Solver::Integrators::EIZ;
  using HIX = typename Solver::Integrators::HIX;
  using HIY = typename Solver::Integrators::HIY;
  using HIZ = typename Solver::Integrators::HIZ;

  static void advance(emdata_t& em) {
    HIX::advance();
    HIY::advance(em.Hy, em.Ez, em.Chyh, em.Ceze);
    HIZ::advance();

    EIX::advance();
    EIY::advance();
    EIZ::advance(em.Ez, em.Hy, em.Jz, em.Ceze, em.Cezh, em.Cjz);
  }
};


struct SolverNone {
  using Integrators = SolverSelector<DefaultSolverArgs,
                                     DefaultSolverArgs,
                                     DefaultSolverArgs,
                                     DefaultSolverArgs,
                                     DefaultSolverArgs,
                                     DefaultSolverArgs>;
};


struct Solver1D {
  using emdata_t = em_data_1d<double>;
  using Integrators = SolverSelector<DefaultSolverArgs,
                                     DefaultSolverArgs,
                                     EzIntegrator_is<Integrator1D<emdata_t, BackwardDiff1D>>,
                                     DefaultSolverArgs,
                                     HyIntegrator_is<Integrator1D<emdata_t, ForwardDiff1D>>,
                                     DefaultSolverArgs>;
};

using Electromagnetics = EM<Solver1D>;

// struct SolverTM {
//   using Integrators = SolverSelector<DefaultSolverArgs,
//                                      DefaultSolverArgs,
//                                      EzIntegrator_is<Integrator2D>,
//                                      HxIntegrator_is<Integrator2D>,
//                                      HyIntegrator_is<Integrator2D>,
//                                      DefaultSolverArgs>;
// };
//
// struct SolverTE {
//   using Integrators = SolverSelector<ExIntegrator_is<Integrator2D>,
//                                      EyIntegrator_is<Integrator2D>,
//                                      DefaultSolverArgs,
//                                      DefaultSolverArgs,
//                                      DefaultSolverArgs,
//                                      HzIntegrator_is<Integrator2D>>;
// };
//
// struct Solver3D {
//   using Integrators = SolverSelector<ExIntegrator_is<Integrator3D>,
//                                      EyIntegrator_is<Integrator3D>,
//                                      EzIntegrator_is<Integrator3D>,
//                                      HxIntegrator_is<Integrator3D>,
//                                      HyIntegrator_is<Integrator3D>,
//                                      HzIntegrator_is<Integrator3D>>;
// };
//
// using SolverTypeList = TypeList<SolverNone, Solver1D, SolverTM, SolverTE, Solver3D>;
//
// using Electromagnetics = EM<TypeListAt<1, SolverTypeList>>;

#endif //EM_SOLVER_H
