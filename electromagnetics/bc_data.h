//
// Created by cepheid on 10/23/24.
//

#ifndef BC_DATA_H
#define BC_DATA_H

// #include "em_traits.h"
// #include "aydenstuff/array.h"
// #include <electromagnetics.h>
// #include "boundaries.h"

#include "electromagnetics.param"

using tf::types::Array1D;
using tf::types::Array2D;
using tf::types::Array3D;

struct NullData {
  explicit NullData(const auto&) {}
};

template<typename Array, EMFace F, EMSide S>
struct PMLData {
  using array_t = Array;
  using value_t = typename Array::value_t;

  explicit PMLData(const Array& f) requires (F == EMFace::X)
  : offsets{get_offsets<F, S, dPML>(f)},
    psi{dPML, f.ny(), f.nz()},
    b{dPML},
    c{dPML}
  {}

  explicit PMLData(const Array& f) requires (F == EMFace::Y)
  : offsets{get_offsets<F, S, dPML>(f)},
    psi{f.nx(), dPML, f.nz()},
    b{dPML},
    c{dPML}
  {}

  explicit PMLData(const Array& f) requires (F == EMFace::Z)
  : offsets{get_offsets<F, S, dPML>(f)},
    psi{f.nx(), f.ny(), dPML},
    b{dPML},
    c{dPML}
  {}

  IntegratorOffsets offsets;
  array_t psi;
  std::vector<value_t> b;
  std::vector<value_t> c;
};

template<typename Array, EMFace F, EMSide S>
struct PeriodicData {
  explicit PeriodicData(const Array& f)
  : offsets{get_offsets<F, S, nHalo>(f)}
  {}

  IntegratorOffsets offsets;
};

template<typename ex_t, typename ey_t, typename ez_t, typename hx_t, typename hy_t, typename hz_t>
struct FaceBCs {
  explicit FaceBCs(const auto& emdata)
  : Ex{emdata.Ex},
    Ey{emdata.Ey},
    Ez{emdata.Ez},
    Hx{emdata.Hx},
    Hy{emdata.Hy},
    Hz{emdata.Hz}
  {}

  ex_t Ex;
  ey_t Ey;
  ez_t Ez;
  hx_t Hx;
  hy_t Hy;
  hz_t Hz;
};

template<typename X0BC, typename X1BC, typename Y0BC, typename Y1BC, typename Z0BC, typename Z1BC>
struct BCData {

  explicit BCData(const auto& emdata)
  : x0{emdata},
    x1{emdata},
    y0{emdata},
    y1{emdata},
    z0{emdata},
    z1{emdata}
  {}

  X0BC x0;
  X1BC x1;
  Y0BC y0;
  Y1BC y1;
  Z0BC z0;
  Z1BC z1;
};

#endif //BC_DATA_H
