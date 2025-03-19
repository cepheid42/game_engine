#ifndef EM_DATA_HPP
#define EM_DATA_HPP

#include "sources.hpp"
#include "bc_data.hpp"
#include "array.hpp"

namespace tf::electromagnetics {
  struct EMData {
    EMData() = delete;
    explicit EMData(std::size_t, std::size_t, std::size_t, double, double);

    void init_coefficients(double, double);

    Array3D<ex_t> Ex;
    Array3D<ex_t> Jx;
    Array3D<ex_t> Cexe;
    Array3D<ex_t> Cexhy;
    Array3D<ex_t> Cexhz;
    Array3D<ex_t> Cjx;

    Array3D<ey_t> Ey;
    Array3D<ey_t> Jy;
    Array3D<ey_t> Ceye;
    Array3D<ey_t> Ceyhx;
    Array3D<ey_t> Ceyhz;
    Array3D<ey_t> Cjy;

    Array3D<ez_t> Ez;
    Array3D<ez_t> Jz;
    Array3D<ez_t> Ceze;
    Array3D<ez_t> Cezhx;
    Array3D<ez_t> Cezhy;
    Array3D<ez_t> Cjz;

    Array3D<hx_t> Hx;
    Array3D<hx_t> Chxh;
    Array3D<hx_t> Chxey;
    Array3D<hx_t> Chxez;

    Array3D<hy_t> Hy;
    Array3D<hy_t> Chyh;
    Array3D<hy_t> Chyex;
    Array3D<hy_t> Chyez;

    Array3D<hz_t> Hz;
    Array3D<hz_t> Chzh;
    Array3D<hz_t> Chzex;
    Array3D<hz_t> Chzey;

    std::vector<CurrentSource> srcs{};

    Array3D<void> empty{}; // for the shits and possibly also some giggles...
  };
} // end namespace tf::electromagnetics



#endif //EM_DATA_HPP
