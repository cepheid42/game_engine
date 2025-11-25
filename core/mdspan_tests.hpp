#ifndef MDSPAN_TESTS_HPP
#define MDSPAN_TESTS_HPP

#include "mdspan.hpp"

#include <vector>
#include <print>
#include <cassert>

constexpr auto nx = 15zu;
constexpr auto ny = 10zu;
constexpr auto nz = 5zu;
constexpr auto size = nx * ny * nz;

// Full 3d span over data
using mdspan_t = std::mdspan<
   std::size_t,
   std::dextents<std::size_t, 3>
>;

   // 3d span
using subspan3d_t = std::mdspan<
   std::size_t,
   std::dextents<std::size_t, 3>,
   std::layout_stride
>;

// 2D span
using subspan2d_t = std::mdspan<
   std::size_t,
   std::dextents<std::size_t, 2>,
   std::layout_stride
>;

// 1D span
using subspan_t = std::mdspan<
   std::size_t,
   std::dextents<std::size_t, 1>,
   std::layout_stride
>;

void fill(auto& md) {
   for (std::size_t i = 0; i < md.extent(0); i++) {
      for (std::size_t j = 0; j < md.extent(1); j++) {
         for (std::size_t k = 0; k < md.extent(2); k++) {
            md[i, j, k] = k  + (nz * j) + (ny * nz * i);
         }
      }
   }
}

void difference(auto& src, auto& dest) {
   for (std::size_t i = 0; i < dest.extent(0); i++) {
      for (std::size_t j = 0; j < dest.extent(1); j++) {
         for (std::size_t k = 0; k < dest.extent(2); k++) {
            dest[i, j, k] = src[i, j, k + 1] - src[i, j, k];
         }
      }
   }
}

void print(auto& md) {
   for (std::size_t i = 0; i < md.extent(0); i++) {
      for (std::size_t j = 0; j < md.extent(1); j++) {
         for (std::size_t k = 0; k < md.extent(2); k++) {
            std::print("{:>3}, ", md[i, j, k]);
         }
         std::println();
      }
      std::println();
   }
   std::println();
}

int test() {
   // Data
   std::vector<std::size_t> test(size);
   std::vector<std::size_t> result_data(size);
   // Using mdspan make indexing subspans easier
   mdspan_t og{test.data(), std::extents{nx, ny, nz}};
   mdspan_t result{result_data.data(), std::extents{nx, ny, nz}};

   fill(og);
   // print(og);

   // /* ---------------------- Backward DX ---------------------- */
   // const subspan3d_t back_dx0{
   //    &og[0, 0, 0],                  // start point
   //    {std::extents{nz, ny, nx},     // slice shape
   //     std::array{1zu, nz, ny * nz}} // slice strides
   // };
   //
   // const subspan3d_t back_dx1{
   //    &result[1, 0, 0],              // start point
   //    {std::extents{nz, ny, nx - 1}, // slice shape
   //     std::array{1zu, nz, ny * nz}} // slice strides
   // };
   // difference(back_dx0, back_dx1);
   // print(result);

   // /* ---------------------- Forward DX ---------------------- */
   // const subspan3d_t forw_dx0{
   //    &og[0, 0, 0],                  // start point
   //    {std::extents{nz, ny, nx},      // slice shape
   //     std::array{1zu, nz, ny * nz}}  // slice strides
   // };
   //
   // const subspan3d_t forw_dx1{
   //    &result[0, 0, 0],              // start point
   //    {std::extents{nz, ny, nx - 1},  // slice shape
   //     std::array{1zu, nz, ny * nz}}  // slice strides
   // };
   // difference(forw_dx0, forw_dx1);
   // print(result);

   // /* ---------------------- Backward DY ---------------------- */
   // const subspan3d_t back_dy0{
   //    &og[0, 0, 0],                  // start point
   //    {std::extents{nx, nz, ny},      // slice shape
   //     std::array{ny * nz, 1zu, nz}}  // slice strides
   // };
   //
   // const subspan3d_t back_dy1{
   //    &result[0, 1, 0],              // start point
   //    {std::extents{nx, nz, ny - 1},  // slice shape
   //     std::array{ny * nz, 1zu, nz}}  // slice strides
   // };
   // difference(back_dy0, back_dy1);
   // print(result);

   // /* ---------------------- Forward DY ---------------------- */
   // const subspan3d_t forw_dy0{
   //    &og[0, 0, 0],                  // start point
   //    {std::extents{nx, nz, ny},      // slice shape
   //     std::array{ny * nz, 1zu, nz}}  // slice strides
   // };
   //
   // const subspan3d_t forw_dy1{
   //    &result[0, 0, 0],              // start point
   //    {std::extents{nx, nz, ny - 1},  // slice shape
   //     std::array{ny * nz, 1zu, nz}}  // slice strides
   // };
   // difference(forw_dy0, forw_dy1);
   // print(result);

   // /* ---------------------- Backward DZ ---------------------- */
   // const subspan3d_t back_dz0{
   //    &og[0, 0, 0],                  // start point
   //    {std::extents{nx, ny, nz},      // slice shape
   //     std::array{ny * nz, nz, 1zu}}  // slice strides
   // };
   //
   // const subspan3d_t back_dz1{
   //    &result[0, 0, 1],              // start point
   //    {std::extents{nx, ny, nz - 1},      // slice shape
   //     std::array{ny * nz, nz, 1zu}}  // slice strides
   // };
   // difference(back_dz0, back_dz1);
   // print(result);

   // /* ---------------------- Forward DZ ---------------------- */
   // const subspan3d_t forw_dz0{
   //    &og[0, 0, 0],                  // start point
   //    {std::extents{nx, ny, nz},      // slice shape
   //     std::array{ny * nz, nz, 1zu}}  // slice strides
   // };
   //
   // const subspan3d_t forw_dz1{
   //    &result[0, 0, 0],              // start point
   //    {std::extents{nx, ny, nz - 1},  // slice shape
   //     std::array{ny * nz, nz, 1zu}}  // slice strides
   // };
   // difference(forw_dz0, forw_dz1);
   // print(result);

   // /* ---------------------- Array Slices ---------------------- */
   // // Constants indices
   // auto x_idx = 2zu;
   // auto y_idx = 2zu;
   // auto z_idx = 2zu;
   //
   // /* ---------------------- XY 2d slice ---------------------- */
   // subspan2d_t xy_span{
   //    &og[0, 0, z_idx],         // Pointer to data start
   //    {std::extents{nx, ny},     // Slice shape
   //     std::array{ny * nz, 1zu}} // Slice strides
   // };
   //
   // for (std::size_t i = 0; i < nx; i++) {
   //    for (std::size_t j = 0; j < ny; j++) {
   //       std::print("{:>3}, ", xy_span[i, j]);
   //    }
   //    std::println();
   // }
   // std::println();
   //
   // /* ---------------------- XZ 2d slice ---------------------- */
   // subspan2d_t xz_span{
   //    &og[0, y_idx, 0],         // Pointer to data start
   //    {std::extents{nx, nz},     // Slice shape
   //     std::array{ny * nz, 1zu}} // Slice strides
   // };
   // for (std::size_t i = 0; i < nx; i++) {
   //    for (std::size_t k = 0; k < nz; k++) {
   //       std::print("{:>3}, ", xz_span[i, k]);
   //    }
   //    std::println();
   // }
   // std::println();
   //
   // /* ---------------------- YZ 2d slice ---------------------- */
   // subspan2d_t yz_span{
   //    &og[x_idx, 0, 0],     // Pointer to data start
   //    {std::extents{ny, nz}, // Slice shape
   //     std::array{nz, 1zu}}  // Slice strides
   // };
   // for (std::size_t j = 0; j < ny; j++) {
   //    for (std::size_t k = 0; k < nz; k++) {
   //       std::print("{:>3}, ", yz_span[j, k]);
   //    }
   //    std::println();
   // }
   // std::println();
   //
   // /* ---------------------- x-axis 1d slice ---------------------- */
   // subspan_t x_slice{
   //    &og[0, y_idx, z_idx], // Pointer to data start
   //    {std::extents{nx},     // Slice shape
   //     std::array{ny * nz}}  // Slice strides
   // };
   //
   // for (std::size_t i = 0; i < nx; i++) {
   //    std::print("{:>3}, ", x_slice[i]);
   // }
   // std::println();
   //
   // /* ---------------------- y-axis 1d slice ---------------------- */
   // subspan_t y_slice{
   //    &og[x_idx, 0, y_idx], // Pointer to data start
   //    {std::extents{ny},     // Slice shape
   //     std::array{nz}}       // Slice strides
   // };
   //
   // for (std::size_t j = 0; j < ny; j++) {
   //    std::print("{:>3}, ", y_slice[j]);
   // }
   // std::println();
   //
   // /* ---------------------- z-axis 1d slice ---------------------- */
   // subspan_t z_slice{
   //    &og[x_idx, y_idx, 0], // Pointer to data start
   //    {std::extents{nz},     // Slice shape
   //     std::array{1zu}}      // Slice strides
   // };
   //
   // for (std::size_t k = 0; k < nz; k++) {
   //    std::print("{:>3}, ", z_slice[k]);
   // }
   // std::println();

   return 0;
}
#endif //MDSPAN_TESTS_HPP
