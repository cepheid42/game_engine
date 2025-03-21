#ifndef GPUEM_TRIDIAGONAL_KERNELS_CUH
#define GPUEM_TRIDIAGONAL_KERNELS_CUH

#include "../Utilities/Utilities.cuh"

namespace TDSKernels {
  __global__ void init_tds_coefs(devMatrix& __restrict__ ac,
                                 devMatrix& __restrict__ b,
                                 const devMatrix& __restrict__ eps_r,
                                 const devMatrix& __restrict__ mu_r,
                                 const devMatrix& __restrict__ sigma_e,
                                 const devMatrix& __restrict__ sigma_m,
                                 const FPTYPE dt,
                                 const FPTYPE da)
  {
    auto i = blockIdx.y;
    auto j = blockIdx.x;
    auto k = threadIdx.x;

    const auto eps = eps_r(i, j, k) * Constants::EPS0;
    const auto mu = mu_r(i, j, k) * Constants::MU0;

    const auto tmp1 = (sigma_e(i, j, k) * dt) / (4.0 * eps);
    const auto tmp2 = (sigma_m(i, j, k) * dt) / (4.0 * mu);
    const auto tmp3 = SQR(dt) / (4.0 * mu * eps * SQR(da));

    const auto alpha = tmp3 / ((1.0 + tmp1) * (1.0 + tmp2));
    const auto beta = 1.0 + (2.0 * alpha);

    ac(i, j, k) = -alpha;
    b(i, j, k) = beta;
  }

  template<FieldComponent C>
  __device__ void load_field_shmem(FPTYPE* __restrict__ d,
                                   const devMatrix& __restrict__ d_d,
                                   const uint32_t systemSize)
  {
    auto bx = blockIdx.y;
    auto by = blockIdx.x;
    auto tid = threadIdx.x;

    #pragma unroll 2
    for (auto i = tid; i < systemSize; i += blockDim.x) {
      if constexpr(C == X) {
        d[i] = d_d(i, bx, by);
      }
      else if constexpr(C == Y) {
        d[i] = d_d(bx, i, by);
      }
      else {
        d[i] = d_d(bx, by, i);
      }
    }

    __syncthreads();
  }

  template<FieldComponent C>
  __device__ void store_field_shmem(devMatrix& __restrict__ d_x,
                                    const FPTYPE* __restrict__ x,
                                    const uint32_t systemSize) {
    auto bx = blockIdx.y;
    auto by = blockIdx.x;
    auto tid = threadIdx.x;

    #pragma unroll 2
    for (auto i = tid; i < systemSize; i += blockDim.x) {
      if constexpr(C == X) {
        d_x(i, bx, by) = x[i];
      }
      else if constexpr(C == Y) {
        d_x(bx, i, by) = x[i];
      }
      else {
        d_x(bx, by, i) = x[i];
      }
    }
  }

  /******** In-Place Cyclic Reduction ********/
  template<FieldComponent C>
  __global__ void CyclicReduction(devMatrix& __restrict__ d_x,
                                  const devMatrix& __restrict__ d_ac,
                                  const devMatrix& __restrict__ d_b)
  {
    extern __shared__ FPTYPE shmem[];

    auto bx = blockIdx.y;
    auto by = blockIdx.x;
    auto tid = threadIdx.x;

    auto nthreads = blockDim.x;
    auto systemSize = d_b.z;

    // Calculate number of divisions by 2 to reduce system
    uint32_t iteration = __float2uint_rn(log2f(__uint2float_rn(nthreads)));

    auto a = (FPTYPE *) shmem;
    auto b = (FPTYPE *) &a[systemSize];
    auto c = (FPTYPE *) &b[systemSize];
    auto d = (FPTYPE *) &c[systemSize];
    auto x = (FPTYPE *) &d[systemSize];

    // Load shared memory
    load_field_shmem<C>(d, d_x, systemSize);

    #pragma unroll 2
    for (uint32_t i = tid; i < systemSize; i += blockDim.x) {
      a[i] = d_ac(bx, by, i);
      c[i] = d_ac(bx, by, i);
      b[i] = d_b(bx, by, i);
      x[i] = 0.0; // shmem is not zero-initialized
    }
    a[0] = 0.0;
    c[systemSize - 1] = 0.0;
    __syncthreads();

    // forward elimination
    int stride = 1;
    int delta;
    for(int j = 0; j < iteration; j++) {
      delta = stride;
      stride *= 2;

      if (tid < nthreads) {
        auto i = (stride * tid) + (stride - 1);
        auto iL = i - delta;
        auto iR = (i + delta >= systemSize) ? systemSize - 1 : i + delta;

        auto k1 = a[i] / b[iL];
        auto k2 = c[i] / b[iR];

        b[i] = b[i] - (k1 * c[iL]) - (k2 * a[iR]);
        d[i] = d[i] - (k1 * d[iL]) - (k2 * d[iR]);
        a[i] = -k1 * a[iL];
        c[i] = -k2 * c[iR];
      }
      nthreads /= 2;
      __syncthreads();
    }

    // Solve remaining two equations
    if (tid == 0) {
      auto a1 = stride - 1;
      auto a2 = 2 * stride - 1;

      auto k3 = (b[a1] * b[a2]) - (c[a1] * a[a2]);

      x[a1] = (b[a2] * d[a1] - c[a1] * d[a2]) / k3;
      x[a2] = (d[a2] * b[a1] - d[a1] * a[a2]) / k3;
    }

    // backward substitution
    nthreads = 2;
    for(int j = 0; j < iteration; j++) {
      delta = stride / 2;

      if (tid < nthreads) {
        auto i = (stride * tid) + (delta - 1);

        if (i == delta - 1) {
          // Keep inbound on left side
          x[i] = (d[i] - (c[i] * x[i + delta])) / b[i];
        } else {
          x[i] = (d[i] - (a[i] * x[i - delta]) - (c[i] * x[i + delta])) / b[i];
        }
      }
      stride = delta;
      nthreads *= 2;
      __syncthreads();
    }

    // Store solution to global memory
    store_field_shmem<C>(d_x, x, systemSize);
  }

  /******** In-Place Parallel Cyclic Reduction ********/
  template<FieldComponent C>
  __global__ void ParallelCyclicReduction(devMatrix& __restrict__ d_x,
                                          const devMatrix& __restrict__ d_ac,
                                          const devMatrix& __restrict__ d_b)
  {
    extern __shared__ FPTYPE shmem[];

    auto bx = blockIdx.y;
    auto by = blockIdx.x;
    auto tid = threadIdx.x;

    auto systemSize = d_b.z;

    // Calculate number of divisions by 2 to reduce system
    auto iteration = __float2uint_rn(log2f(__uint2float_rn(blockDim.x))) - 1;

    auto a = (FPTYPE *) shmem;
    auto b = (FPTYPE *) &a[systemSize];
    auto c = (FPTYPE *) &b[systemSize];
    auto d = (FPTYPE *) &c[systemSize];
    auto x = (FPTYPE *) &d[systemSize];

    // Load shared memory
    load_field_shmem<C>(d, d_x, systemSize);

    a[tid] = d_ac(bx, by, tid);
    c[tid] = d_ac(bx, by, tid);
    b[tid] = d_b(bx, by, tid);
    x[tid] = 0.0; // shmem is not zero-initialized

    a[0] = 0.0;
    c[systemSize - 1] = 0.0;
    __syncthreads();

    FPTYPE aNew, bNew, cNew, dNew;
    auto delta = 1;
    for (int j = 0; j < iteration; j++) {
        int i = tid;

        int iR = i + delta;
        if (iR >= systemSize) {
            iR = systemSize - 1;
        }

        int iL= i - delta;
        if (iL < 0) {
            iL = 0;
        }

        auto k1 = a[i] / b[iL];
        auto k2 = c[i] / b[iR];

        bNew = b[i] - (k1 * c[iL]) - (k2 * a[iR]);
        dNew = d[i] - (k1 * d[iL]) - (k2 * d[iR]);
        aNew = -k1 * a[iL];
        cNew = -k2 * c[iR];
        __syncthreads();

        a[i] = aNew;
        b[i] = bNew;
        c[i] = cNew;
        d[i] = dNew;
        __syncthreads();

        delta *= 2;
    }

    if (tid < delta) {
        auto a1 = tid;
        auto a2 = tid + delta;

        auto tmp3 = (b[a2] * b[a1]) - (c[a1] * a[a2]);
        x[a1] = ((b[a2] * d[a1]) - (c[a1] * d[a2])) / tmp3;
        x[a2] = ((d[a2] * b[a1]) - (d[a1] * a[a2])) / tmp3;
    }
    __syncthreads();

    // Store solution to global memory
    store_field_shmem<C>(d_x, x, systemSize);
  }
}
#endif //GPUEM_TRIDIAGONAL_KERNELS_CUH
