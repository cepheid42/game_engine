//
// Created by cepheid on 9/23/24.
//

#ifndef EM_SOLVER_H
#define EM_SOLVER_H

#include <cassert>
#include <iostream>

#include "./aydenstuff/tags.h"
// #include "em_data.h"


// template<typename EMClass>
// struct EIntegrator1D {
//   static void advanceE() {
//     std::cout << "Hi, EIntegrator1D." << std::endl;
//     // EMClass::advanceE1D("EIntegrator");
//   }
// };
//
// template<typename EMClass>
// struct EIntegrator2D {
//   static void advanceE() {
//     std::cout << "Hi, EIntegrator2D." << std::endl;
//     // EMClass::advanceE1D("EIntegrator");
//   }
// };
//
// template<typename EMClass>
// struct BIntegrator1D {
//   static void advanceB() {
//     std::cout << "Hi, BIntegrator1D." << std::endl;
//     // EMClass::advanceB1D("BIntegrator");
//   }
// };
//
// template<typename EMClass>
// struct BIntegrator2D {
//   static void advanceB() {
//     std::cout << "Hi, BIntegrator2D." << std::endl;
//     // EMClass::advanceB1D("BIntegrator");
//   }
// };

enum class Fields { Ex, Ey, Ez, Hx, Hy, Hz };

template<Fields Field>
struct ExplicitFunctor1D {
  static void apply(size_t i) {
    if constexpr (Field == Fields::Ez) {
      std::cout << "ExplicitFunctor::Ez(" << i <<")" << std::endl;
    } else {
      std::cout << "ExplicitFunctor::Hy(" << i <<")" << std::endl;
    }
  }
};


template<typename EMClass, typename Func>
struct Integrator1D {
  static void advance() {
    for (size_t i = 0; i < 2; i++) {
      Func::apply(i);
    }
  }
};

using ElectricUpdateEz = ExplicitFunctor1D<Fields::Ez>;
using MagneticUpdateHy = ExplicitFunctor1D<Fields::Hy>;

struct NullIntegrator {
  static void advance() { std::cout << "NullIntegrator::advance()" << std::endl; }
};

template<typename EMClass>
struct Solver1D : Integrator1D<EMClass, ElectricUpdateEz>, Integrator1D<EMClass, MagneticUpdateHy> {
  using EIX = NullIntegrator;
  using EIY = NullIntegrator;
  using EIZ = Integrator1D<EMClass, ElectricUpdateEz>;
  using HIX = NullIntegrator;
  using HIY = Integrator1D<EMClass, MagneticUpdateHy>;
  using HIZ = NullIntegrator;
};

// template<typename EMClass>
// struct Solver2D : EIntegrator2D<EMClass>, BIntegrator2D<EMClass> {};




template<template<typename> typename Solver>
struct EM : Solver<EM<Solver>> {
  using EIX = typename Solver<EM>::EIX;
  using EIY = typename Solver<EM>::EIY;
  using EIZ = typename Solver<EM>::EIZ;
  using HIX = typename Solver<EM>::HIX;
  using HIY = typename Solver<EM>::HIY;
  using HIZ = typename Solver<EM>::HIZ;

  static void advance() {
    HIX::advance();
    HIY::advance();
    HIZ::advance();

    EIX::advance();
    EIY::advance();
    EIZ::advance();
  }
};

using EM1D = EM<Solver1D>;
// using EM2D = EM<Solver2D>;


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
