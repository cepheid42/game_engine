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

//============== Integrator Classes =============
//===============================================
struct Integrator1D {
  static void advance() { std::cout << "Integrator1D::advance()" << std::endl; }
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
  using EIX = typename Solver::Integrators::EIX;
  using EIY = typename Solver::Integrators::EIY;
  using EIZ = typename Solver::Integrators::EIZ;
  using HIX = typename Solver::Integrators::HIX;
  using HIY = typename Solver::Integrators::HIY;
  using HIZ = typename Solver::Integrators::HIZ;

  static void advance() {
    HIX::advance();
    HIY::advance();
    HIZ::advance();

    EIX::advance();
    EIY::advance();
    EIZ::advance();
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
  using Integrators = SolverSelector<DefaultSolverArgs,
                                     DefaultSolverArgs,
                                     EzIntegrator_is<Integrator1D>,
                                     DefaultSolverArgs,
                                     HyIntegrator_is<Integrator1D>,
                                     DefaultSolverArgs>;
};

struct SolverTM {
  using Integrators = SolverSelector<DefaultSolverArgs,
                                     DefaultSolverArgs,
                                     EzIntegrator_is<Integrator2D>,
                                     HxIntegrator_is<Integrator2D>,
                                     HyIntegrator_is<Integrator2D>,
                                     DefaultSolverArgs>;
};

struct SolverTE {
  using Integrators = SolverSelector<ExIntegrator_is<Integrator2D>,
                                     EyIntegrator_is<Integrator2D>,
                                     DefaultSolverArgs,
                                     DefaultSolverArgs,
                                     DefaultSolverArgs,
                                     HzIntegrator_is<Integrator2D>>;
};

struct Solver3D {
  using Integrators = SolverSelector<ExIntegrator_is<Integrator3D>,
                                     EyIntegrator_is<Integrator3D>,
                                     EzIntegrator_is<Integrator3D>,
                                     HxIntegrator_is<Integrator3D>,
                                     HyIntegrator_is<Integrator3D>,
                                     HzIntegrator_is<Integrator3D>>;
};

using SolverTypeList = TypeList<SolverNone, Solver1D, SolverTM, SolverTE, Solver3D>;

using Electromagnetics = EM<TypeListAt<3, SolverTypeList>>;

#endif //EM_SOLVER_H
