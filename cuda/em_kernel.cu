#include "em_kernel.cuh"

#include <iostream>


template struct cudaClass<Thingy>;

#define cudaChk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char* file, int line, bool abort=true) {
   if (code != cudaSuccess) {
      std::cout << "[" << file << ":" << line << "] GPU Error: " << cudaGetErrorString(code) << std::endl;
      if (abort) { exit(code); }
   }
}

template<typename T>
void cudaClass<T>::run() const {
   float* result;
   cudaChk(cudaMallocManaged(&result, 32 * sizeof(float)));
   cudaChk(cudaDeviceSynchronize());


   kernel<T><<<1, 32>>>(t, 2.0f, result);
   cudaChk(cudaDeviceSynchronize());

   for (int i = 0; i < 32; i++) {
      std::cout << result[i] << ", ";
   }
   std::cout << std::endl;

   cudaChk(cudaFree(result));
   cudaChk(cudaDeviceSynchronize());
}


