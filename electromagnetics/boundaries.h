//
// Created by cepheid on 10/17/24.
//

#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include <iostream>

//====== 1D Boundaries ========
//=============================
template<typename T>
struct Boundary {
  void apply() {
    static_cast<T*>(this)->apply();
  }
};

template<typename Array>
struct Periodic1D : Boundary<Periodic1D<Array>> {
  using value_t = typename Array::value_t;

  void apply() {
    std::cout << "Periodic1D" << std::endl;
    // f[0] = f[f.nx - 1];
    // f[f.nx - 2] = f[1];
  }
};

// template<typename BC>
// struct PeriodicBoundary {
//
// };



//====== 2D Boundaries ========
//=============================


//====== 3D Boundaries ========
//=============================

#endif //BOUNDARIES_H
