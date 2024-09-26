//
// Created by cepheid on 9/23/24.
//

#ifndef EM_SOLVER_H
#define EM_SOLVER_H

#include <cassert>
#include <iostream>

#include "./aydenstuff/tags.h"
#include "em_data.h"


enum class Fields { Ex, Ey, Ez, Hx, Hy, Hz };

template<template<typename> typename EMData, Fields Field>
struct ExplicitUpdate;

template<Fields Field>
struct ExplicitUpdate<em_data_1d, Field> {
  using emdata_t = em_data_1d<tf::types::Array1D<double>>;

  static void apply(emdata_t& em, const size_t i) {
    if constexpr (Field == Fields::Ez) {
      std::cout << "ExplicitFunctor1D::Ez(" << i << ") = " << em.Ez.size() << std::endl;
    } else {
      std::cout << "ExplicitFunctor1D::Hy(" << i << ") = " << em.Hy.size() << std::endl;
    }
  }
};

template<Fields Field>
struct ExplicitUpdate<em_data_tm, Field> {
  using emdata_t = em_data_tm<tf::types::Array2D<double>>;

  static void apply(emdata_t& em, size_t i, size_t j) {
    if constexpr (Field == Fields::Ez) {
      std::cout << "ExplicitFunctorTM_2D::Ez(" << i << ", " << j << ") = " << em.Ez.nx << ", " << em.Ez.nz << std::endl;
    } else if constexpr (Field == Fields::Hx) {
      std::cout << "ExplicitFunctor2D::TM_Hx(" << i << ", " << j << ") = " << em.Hx.nx << ", " << em.Hx.nz << std::endl;
    } else {
      std::cout << "ExplicitFunctor2D::TM_Hy(" << i << ", " << j << ") = " << em.Hx.nx << ", " << em.Hx.nz << std::endl;
    }
  }
};


template<typename EMClass, typename Func>
struct Integrator1D {
  using emdata_t = em_data_1d<tf::types::Array1D<double>>;

  static void advance(emdata_t& em) {
    for (size_t i = 0; i < 2; i++) {
      Func::apply(em, i);
    }
  }
};

template<typename EMClass, typename Func>
struct Integrator2D {
  using emdata_t = em_data_tm<tf::types::Array2D<double>>;

  static void advance(emdata_t& em) {
    for (size_t i = 0; i < 2; i++) {
      for (size_t j = 0; j < 2; j++) {
        Func::apply(em, i, j);
      }
    }
  }
};

template<template<typename> typename EMData>
using ElectricUpdateEx = ExplicitUpdate<EMData, Fields::Ex>;
template<template<typename> typename EMData>
using ElectricUpdateEy = ExplicitUpdate<EMData, Fields::Ey>;
template<template<typename> typename EMData>
using ElectricUpdateEz = ExplicitUpdate<EMData, Fields::Ez>;
template<template<typename> typename EMData>
using MagneticUpdateHx = ExplicitUpdate<EMData, Fields::Hx>;
template<template<typename> typename EMData>
using MagneticUpdateHy = ExplicitUpdate<EMData, Fields::Hy>;
template<template<typename> typename EMData>
using MagneticUpdateHz = ExplicitUpdate<EMData, Fields::Hz>;

struct NullIntegrator {
  static void advance(auto&) { std::cout << "NullIntegrator::advance()" << std::endl; }
};


template<typename EMClass>
struct Solver1D
: Integrator1D<EMClass, ElectricUpdateEz<em_data_1d>>,
  Integrator1D<EMClass, MagneticUpdateHy<em_data_1d>>
{
  using emdata_t = em_data_1d<tf::types::Array1D<double>>;

  using EIX = NullIntegrator;
  using EIY = NullIntegrator;
  using EIZ = Integrator1D<EMClass, ElectricUpdateEz<em_data_1d>>;
  using HIX = NullIntegrator;
  using HIY = Integrator1D<EMClass, MagneticUpdateHy<em_data_1d>>;
  using HIZ = NullIntegrator;
};

template<typename EMClass>
struct SolverTM
: Integrator2D<EMClass, ElectricUpdateEz<em_data_tm>>,
  Integrator2D<EMClass, MagneticUpdateHx<em_data_tm>>,
  Integrator2D<EMClass, MagneticUpdateHy<em_data_tm>>
{
  using emdata_t = em_data_tm<tf::types::Array2D<double>>;

  using EIX = NullIntegrator;
  using EIY = NullIntegrator;
  using EIZ = Integrator2D<EMClass, ElectricUpdateEz<em_data_tm>>;
  using HIX = Integrator2D<EMClass, MagneticUpdateHx<em_data_tm>>;
  using HIY = Integrator2D<EMClass, MagneticUpdateHy<em_data_tm>>;
  using HIZ = NullIntegrator;
};



template<template<typename> typename Solver>
struct EM : Solver<EM<Solver>> {
  using EIX = typename Solver<EM>::EIX;
  using EIY = typename Solver<EM>::EIY;
  using EIZ = typename Solver<EM>::EIZ;
  using HIX = typename Solver<EM>::HIX;
  using HIY = typename Solver<EM>::HIY;
  using HIZ = typename Solver<EM>::HIZ;

  using emdata_t = typename Solver<EM>::emdata_t;

  static void advance(emdata_t& em) {
    HIX::advance(em);
    HIY::advance(em);
    HIZ::advance(em);

    EIX::advance(em);
    EIY::advance(em);
    EIZ::advance(em);
  }
};

using EM1D = EM<Solver1D>;
using EMTM = EM<SolverTM>;


// using X12 = X<ExtraFeature1, ExtraFeature2>;

// //============== Solver Functors =============
// //============================================
// // template<typename T>
// // struct explicit_solver {
// //   static void apply()
// // };
//
//
// //========= Field Integrator Classes =========
// //============================================
// template<typename Derived>
// struct ElectricIntegrator;
//
// // template<>
// // struct ElectricIntegrator<em_data_1d> {
// //   // using array_t = tf::types::Array1D<T>;
// //   //
// //   // static void apply(array_t& f, const array_t& d1, const array_t& d2, const array_t& j, const array_t& c_f, const array_t& c_d, const array_t& c_j) {
// //   //   for (std::size_t i = 0; i < f.nx; i++) {
// //   //
// //   //   }
// //   // }
// // };
//
//
//
//
// //====== Solver Specialization Classes =======
// //============================================
// template<typename EMClass>
// struct Solver_1D : ElectricIntegrator<EMClass> {
//   void advance() {}
// };
//
// // template<template<typename> typename EMClass>
// // struct Solver_TM : Integrator<EMClass>, Integrator<EMClass> {};
// //
// // template<template<typename> typename EMClass>
// // struct Solver_TE : Integrator<EMClass>, Integrator<EMClass> {};
// //
// // template<template<typename> typename EMClass>
// // struct Solver_3D : Integrator<EMClass>, Integrator<EMClass> {};
//
//
//
// //====== Electromagnetic Superclass =======
// //=========================================
// template<typename Solver>
// struct Electromagnetics : Solver<Electromagnetics<Solver>> {
//
//
// };
//
// // template<typename T>
// // using Electromagnetics1D = Electromagnetics<T, Solver_1D<T, em_data_1d>>;
//
// // using ElectromagneticsTM = Electromagnetics<Solver_TM<em_data_tm>>;
// // using ElectromagneticsTE = Electromagnetics<Solver_TE<em_data_te>>;
// // using Electromagnetics3D = Electromagnetics<Solver_3D<em_data_3d>>;
#endif //EM_SOLVER_H
