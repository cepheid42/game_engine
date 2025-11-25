#ifndef EM_ARRAY_HPP
#define EM_ARRAY_HPP

#include "vec3.hpp"
#include "traits.hpp"

#include <vector>

namespace tf {
template<typename T>
class Array3D {
   using value_t = T;

public:
   explicit Array3D() = default;

   explicit Array3D(const std::size_t nx, const std::size_t ny, const std::size_t nz)
   : nx_(nx),
     ny_(ny),
     nz_(nz),
     data_(nx * ny * nz)
   {}

   explicit Array3D(const vec3<T>& dims_) : Array3D(dims_.x, dims_.y, dims_.z) {}

   auto begin() { return data_.begin(); }
   auto   end() { return data_.end(); }
   auto* data() { return data_.data(); }
   
   [[nodiscard]] auto      size() const -> std::size_t { return data_.size(); }
   [[nodiscard]] auto num_bytes() const -> std::size_t { return size() * sizeof(T); }
   [[nodiscard]] auto get_scid(const auto i, const auto j, const auto k) const -> std::size_t { return k + nz_ * j + nz_ * ny_ * i; }

   // Specialized accessors
   auto operator()(const auto i, const auto j, const auto k)       ->       T& { return data_[get_scid(i, j, k)]; }
   auto operator()(const auto i, const auto j, const auto k) const -> const T& { return data_[get_scid(i, j, k)]; }

   auto operator[](const auto i)       ->       T& { return data_[i]; }
   auto operator[](const auto i) const -> const T& { return data_[i]; }

   void fill(T value) { for (auto& el: data_) el = value; }

   // [[nodiscard]] bool is_inbounds(const std::size_t i, const std::size_t j, const std::size_t k) const { return i < nx_ and j < ny_ and k < nz_; }

   // Dims
   [[nodiscard]] auto dims() const -> vec3<std::size_t> { return {nx_, ny_, nz_}; }
   [[nodiscard]] auto nx() const -> std::size_t { return nx_; }
   [[nodiscard]] auto ny() const -> std::size_t { return ny_; }
   [[nodiscard]] auto nz() const -> std::size_t { return nz_; }

private:
   // Stride data
   std::size_t    nx_{};
   std::size_t    ny_{};
   std::size_t    nz_{};
   std::vector<T> data_{};
}; // end class Array3D

struct null_t {};

template<>
class Array3D<null_t> {
   using value_t = null_t;
public:
   explicit Array3D() = default;
   explicit Array3D(const std::size_t, const std::size_t, const std::size_t, const auto) {}
   // Specialized accessors
   constexpr auto operator()(const std::size_t, const std::size_t, const std::size_t) const { return 0.0f; }
   constexpr auto operator[](const std::size_t) const { return 0.0f; }
   static constexpr void fill(auto) {}
};
} // end namespace tf

#endif //EM_ARRAY_HPP
