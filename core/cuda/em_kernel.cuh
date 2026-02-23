#ifndef EM_KERNEL_CUH
#define EM_KERNEL_CUH

#include <cuda_runtime.h>
#include <device_launch_parameters.h>

struct Thingy {
   __device__ static float apply(const float x) { return x * x; }
};

template<typename T>
__global__ void kernel(const T t, const float value, float* results) {
   const auto tid = threadIdx.x;
   results[tid] = t.apply(value);
}

template<typename T>
struct cudaClass {
   const T& t;

   explicit cudaClass(const T& t) : t(t) {}

   void run() const;
};

#endif //EM_KERNEL_CUH
