//
// Created by cepheid on 9/23/24.
//

#ifndef EM_SOLVER_H
#define EM_SOLVER_H

#include "em_data.h"

struct TMFunctor {};

template<template<typename> typename EMClass, int Field>
struct Integrator;

template<>
struct Integrator<em_data_tm, 0> {};

template<>
struct Integrator<em_data_tm, 1> {};

template<template<typename> typename EMClass>
struct Solver_TM : Integrator<EMClass, 0>, Integrator<EMClass, 1> {};

template<template<typename> typename EMClass>
struct Solver_TE : Integrator<EMClass, 2>, Integrator<EMClass, 3> {};


template<typename Solver>
struct Electromagnetics : Solver<Electromagnetics<Solver>> {

};

using ElectromagneticsTM = Electromagnetics<Solver_TM<em_data_tm>>;
using ElectromagneticsTE = Electromagnetics<Solver_TE<em_data_te>>;

#endif //EM_SOLVER_H
