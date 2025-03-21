#ifndef EM_DATA_HPP
#define EM_DATA_HPP

#include "program_params.hpp"
#include "sources.hpp"
#include "array.hpp"

namespace tf::electromagnetics {
  struct EMData {

    EMData() = delete;
    explicit EMData(std::size_t, std::size_t, std::size_t, compute_t, compute_t);

    void init_coefficients(compute_t, compute_t);

    Array3D<compute_t> Ex;
    Array3D<compute_t> Jx;
    Array3D<compute_t> Cexe;
    Array3D<compute_t> Cexhy;
    Array3D<compute_t> Cexhz;
    Array3D<compute_t> Cjx;

    Array3D<compute_t> Ey;
    Array3D<compute_t> Jy;
    Array3D<compute_t> Ceye;
    Array3D<compute_t> Ceyhx;
    Array3D<compute_t> Ceyhz;
    Array3D<compute_t> Cjy;

    Array3D<compute_t> Ez;
    Array3D<compute_t> Jz;
    Array3D<compute_t> Ceze;
    Array3D<compute_t> Cezhx;
    Array3D<compute_t> Cezhy;
    Array3D<compute_t> Cjz;

    Array3D<compute_t> Hx;
    Array3D<compute_t> Chxh;
    Array3D<compute_t> Chxey;
    Array3D<compute_t> Chxez;

    Array3D<compute_t> Hy;
    Array3D<compute_t> Chyh;
    Array3D<compute_t> Chyex;
    Array3D<compute_t> Chyez;

    Array3D<compute_t> Hz;
    Array3D<compute_t> Chzh;
    Array3D<compute_t> Chzex;
    Array3D<compute_t> Chzey;

    std::vector<CurrentSource> srcs{};

    Array3D<void> empty{}; // for the shits and possibly also some giggles...
  };
} // end namespace tf::electromagnetics



#endif //EM_DATA_HPP
