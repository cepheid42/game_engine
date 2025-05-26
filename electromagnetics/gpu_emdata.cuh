#ifndef GPU_EMDATA_CUH
#define GPU_EMDATA_CUH

#include <cuda_runtime.h>

namespace tf::gpu
{

template<typename T>
void init_data_pointers(auto& emdata, const std::size_t nx, const std::size_t ny, const std::size_t nz)
{
   cudaMallocManaged(emdata.Ex, sizeof(T) * (nx - 1) * ny * nz);
   cudaMallocManaged(emdata.Jx, sizeof(T) * (nx - 1) * ny * nz);
   cudaMallocManaged(emdata.Cexe, sizeof(T) * (nx - 1) * ny * nz);
   cudaMallocManaged(emdata.Cexhy, sizeof(T) * (nx - 1), ny, nz);
   cudaMallocManaged(emdata.Cexhz, sizeof(T) * (nx - 1), ny, nz);
   cudaMallocManaged(emdata.Cjx, sizeof(T) * (nx - 1) * ny * nz);
   cudaMallocManaged(emdata.Ey, sizeof(T) * nx, (ny - 1), nz));
   cudaMallocManaged(emdata.Jy, sizeof(T) * nx, (ny - 1), nz));
   cudaMallocManaged(emdata.Ceye, sizeof(T) * nx, (ny - 1), nz)),
   cudaMallocManaged(emdata.Ceyhx, sizeof(T) * nx, (ny - 1), nz)),
   cudaMallocManaged(emdata.Ceyhz, sizeof(T) * nx, (ny - 1), nz)),
   cudaMallocManaged(emdata.Cjy, sizeof(T) * nx, (ny - 1), nz)),
   cudaMallocManaged(emdata.Ez, sizeof(T) * nx, ny, (nz - 1)),
   cudaMallocManaged(emdata.Jz, sizeof(T) * nx, ny, (nz - 1)),
   cudaMallocManaged(emdata.Ceze, sizeof(T) * nx, ny, nz - 1)),
   cudaMallocManaged(emdata.Cezhx, sizeof(T) * nx, ny, nz - 1)),
   cudaMallocManaged(emdata.Cezhy, sizeof(T) * nx, ny, nz - 1)),
   cudaMallocManaged(emdata.Cjz, sizeof(T) * nx, ny, nz - 1)),
   cudaMallocManaged(emdata.Hx, sizeof(T) * nx, ny - 1, nz - 1)),
   cudaMallocManaged(emdata.Chxh, sizeof(T) * nx, ny - 1, nz - 1)),
   cudaMallocManaged(emdata.Chxey, sizeof(T) * nx, ny - 1, nz - 1)),
   cudaMallocManaged(emdata.Chxez, sizeof(T) * nx, ny - 1, nz - 1)),
   cudaMallocManaged(emdata.Hy, sizeof(T) * (nx - 1, ny, nz - 1)),
   cudaMallocManaged(emdata.Chyh, sizeof(T) * (nx - 1, ny, nz - 1)),
   cudaMallocManaged(emdata.Chyex, sizeof(T) * (nx - 1, ny, nz - 1)),
   cudaMallocManaged(emdata.Chyez, sizeof(T) * (nx - 1, ny, nz - 1)),
   cudaMallocManaged(emdata.Hz, sizeof(T) * (nx - 1, ny - 1, nz)),
   cudaMallocManaged(emdata.Chzh, sizeof(T) * (nx - 1, ny - 1, nz)),
   cudaMallocManaged(emdata.Chzex, sizeof(T) * (nx - 1, ny - 1, nz)),
   cudaMallocManaged(emdata.Chzey, sizeof(T) * (nx - 1, ny - 1, nz)),
   cudaMallocManaged(emdata.Bx, sizeof(T) * (nx, ny - 1, nz - 1)),
   cudaMallocManaged(emdata.Chxey2, sizeof(T) * (nx, ny - 1, nz - 1)),
   cudaMallocManaged(emdata.Chxez2, sizeof(T) * (nx, ny - 1, nz - 1)),
   cudaMallocManaged(emdata.By, sizeof(T) * (nx - 1, ny, nz - 1)),
   cudaMallocManaged(emdata.Chyex2, sizeof(T) * (nx - 1, ny, nz - 1)),
   cudaMallocManaged(emdata.Chyez2, sizeof(T) * (nx - 1, ny, nz - 1)),
   cudaMallocManaged(emdata.Bz, sizeof(T) * (nx - 1, ny - 1, nz)),
   cudaMallocManaged(emdata.Chzex2, sizeof(T) * (nx - 1, ny - 1, nz)),
   cudaMallocManaged(emdata.Chzey2, sizeof(T) * (nx - 1, ny - 1, nz)),
   cudaMallocManaged(emdata.Ex_app, sizeof(T) * (nx - 1, ny, nz)),
   cudaMallocManaged(emdata.Ey_app, sizeof(T) * (nx, ny - 1, nz)),
   cudaMallocManaged(emdata.Ez_app, sizeof(T) * (nx, ny, nz - 1)),
   cudaMallocManaged(emdata.Bx_app, sizeof(T) * (nx, ny - 1, nz - 1)),
   cudaMallocManaged(emdata.By_app, sizeof(T) * (nx - 1, ny, nz - 1)),
   cudaMallocManaged(emdata.Bz_app, sizeof(T) * (nx - 1, ny - 1, nz))
}

template<typename T>
struct EMData
{
   EMData(const std::size_t nx, const std::size_t ny, const std::size_t nz, const T cfl, const T dt)

   {
      init_coefficients(cfl, dt);
   }

   T* Ex;
   T* Jx;
   T* Cexe;
   T* Cexhy;
   T* Cexhz;
   T* Cjx;

   T* Ey;
   T* Jy;
   T* Ceye;
   T* Ceyhx;
   T* Ceyhz;
   T* Cjy;

   T* Ez;
   T* Jz;
   T* Ceze;
   T* Cezhx;
   T* Cezhy;
   T* Cjz;

   T* Hx;
   T* Chxh;
   T* Chxey;
   T* Chxez;

   T* Hy;
   T* Chyh;
   T* Chyex;
   T* Chyez;

   T* Hz;
   T* Chzh;
   T* Chzex;
   T* Chzey;

   T* Bx;
   T* Chxey2;
   T* Chxez2;

   T* By;
   T* Chyex2;
   T* Chyez2;

   T* Bz;
   T* Chzex2;
   T* Chzey2;

   T* Ex_app;
   T* Ey_app;
   T* Ez_app;

   T* Bx_app;
   T* By_app;
   T* Bz_app;
};

}


#endif //GPU_EMDATA_CUH
