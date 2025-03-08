#ifndef EM_DATA_HPP
#define EM_DATA_HPP

#include <cstdint>

#include "em_array.hpp"

// todo: move these elsewhere
constexpr std::size_t Nx = 100lu;
constexpr std::size_t Ny = 100lu;
constexpr std::size_t Nz = 100lu;

namespace tf::electromagnetics {

  class EMData {
    using array_t = Array3D<double, Nx, Ny, Nz>;
    
    EMData() = delete;
    explicit EMData(int);
    void init_coefficients();

    array_t Ex;
    array_t Jx;
    array_t Cexe;
    array_t Cexhy;
    array_t Cexhz;
    array_t Cjx;

    array_t Ey;
    array_t Jy;
    array_t Ceye;
    array_t Ceyhx;
    array_t Ceyhz;
    array_t Cjy;

    array_t Ez;
    array_t Jz;
    array_t Ceze;
    array_t Cezhx;
    array_t Cezhy;
    array_t Cjz;

    array_t Hx;
    array_t Chxey;
    array_t Chxez;
    array_t Chxh;

    array_t Hy;
    array_t Chyex;
    array_t Chyez;
    array_t Chyh;

    array_t Hz;
    array_t Chzex;
    array_t Chzey;
    array_t Chzh;
  };
}

#endif //EM_DATA_HPP
