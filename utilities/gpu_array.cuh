#ifndef EM_ARRAY_HPP
#define EM_ARRAY_HPP


#include <cuda_runtime.h>

namespace tf::gpu
{
struct Array3D {
   dim3 dims;


   // Constructor
   Array3D(const std::size_t nx_, const std::size_t ny, const std::size_t nz_) {
      dims = make_uint3(nx_, ny_, nz_);
      cudaChk(cudaMallocManaged(&data, Nx * Ny * Nz * sizeof(float)))
    }

   // Copy Constructor
   Array3D(const Array3D &m) {
      dims = make_uint3(m.dims.x, m.dims.y, m.dims.z);
      cudaChk(cudaMallocManaged(&data, m.get_size() * sizeof(float)))
      memcpy(data, m.data, m.get_size() * sizeof(float));
   }

   __device__ __host__ inline uint get_size() const {
      return dims.x * dims.y * dims.z;
   };


   __device__ __host__ inline uint get_index(uint i, uint j, uint k) const {
      return (dims.z * dims.y * i) + (dims.z * j) + k;
   }

   void prefetch_to_host() const {
      cudaChk(cudaMemPrefetchAsync(data, get_size() * sizeof(float), cudaCpuDeviceId))
    }

   void prefetch_to_device() const {
      cudaChk(cudaMemPrefetchAsync(data, get_size() * sizeof(float), 0))
    }

   void memcpyasync_to_host(float *buffer) const {
      cudaChk(cudaMemcpyAsync(buffer, data, get_size() * sizeof(float), cudaMemcpyDeviceToHost))
    }

   void memcpy_to_host(float *buffer) const {
      cudaChk(cudaMemcpy(buffer, data, get_size() * sizeof(float), cudaMemcpyDeviceToHost))
    }
};
} // end namespace tf

#endif //EM_ARRAY_HPP
