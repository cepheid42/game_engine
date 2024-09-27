//
// Created by cepheid on 7/4/24.
//

#ifndef FIELDS_H
#define FIELDS_H

#include <vector>
using fptype = double;
using fp_arr = std::vector<double>;

//=================== Field Classes ========================
//==========================================================
struct Fields1D
{
    explicit Fields1D(const size_t nx)
    : Ez(nx, 0.0),
      Jz(nx, 0.0),
      Hy(nx, 0.0)
    {}

    fp_arr Ez;
    fp_arr Jz;
    fp_arr Hy;
};

struct Fields2D
{
    explicit Fields2D(const size_t nx, const size_t ny)
    : Ez(nx * ny, 0.0),
      Jz(nx * ny, 0.0),
      Hx(nx * (ny - 1), 0.0),
      Hy((nx - 1) * ny, 0.0)
    {}

    fp_arr Ez;
    fp_arr Jz;
    fp_arr Hx;
    fp_arr Hy;
};

struct Fields3D
{
    explicit Fields3D(const size_t nx, const size_t ny, const size_t nz)
    : Ex((nx - 1) * ny * nz, 0.0),
      Ey(nx * (ny - 1) * nz, 0.0),
      Ez(nx * ny * (nz - 1), 0.0),
      Jx((nx - 1) * ny * nz, 0.0),
      Jy(nx * (ny - 1) * nz, 0.0),
      Jz(nx * ny * (nz - 1), 0.0),
      Hx(nx * (ny - 1) * (nz - 1), 0.0),
      Hy((nx - 1) * ny * (nz - 1), 0.0),
      Hz((nx - 1) * (ny - 1) * nz, 0.0)
    {}

    fp_arr Ex;
    fp_arr Ey;
    fp_arr Ez;
    fp_arr Jx;
    fp_arr Jy;
    fp_arr Jz;
    fp_arr Hx;
    fp_arr Hy;
    fp_arr Hz;
};



#endif //FIELDS_H
