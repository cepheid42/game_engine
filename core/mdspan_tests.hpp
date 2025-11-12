#ifndef MDSPAN_TESTS_HPP
#define MDSPAN_TESTS_HPP

#include "mdspan.hpp"

#include <vector>
#include <print>

inline void showcase_mdspan() {
   std::size_t nx = 5uz;
   std::size_t ny = 5uz;
   std::size_t nz = 5uz;

   // Full 3d span over data
   using mdspan_t = std::mdspan<
      int,
      std::dextents<std::size_t, 3>
   >;

   // 2D span
   using subspan2d_t = std::mdspan<
      int,
      std::dextents<std::size_t, 2>,
      std::layout_stride
   >;

   // 1D span
   using subspan_t = std::mdspan<
      int,
      std::dextents<std::size_t, 1>,
      std::layout_stride
   >;

   // Data
   std::vector<int> test(nx*ny*nz);
   for (int i = 0; i < nx*ny*nz; ++i) {
      test[i] = i;
   }

   // Using mdspan make indexing subspans easier
   mdspan_t md1{test.data(), std::extents{nx, ny, nz}};

   // Using mdspan to fill array
   // for (std::size_t i = 0; i < nx; i++) {
   //    for (std::size_t j = 0; j < ny; j++) {
   //       for (std::size_t k = 0; k < nz; k++) {
   //          md1[i, j, k] = k;
   //       }
   //    }
   // }

   // Printing
   // for (std::size_t i = 0; i < nx; i++) {
   //    for (std::size_t j = 0; j < ny; j++) {
   //       for (std::size_t k = 0; k < nz; k++) {
   //          std::print("{:>3}, ", md1[i, j, k]);
   //       }
   //       std::println();
   //    }
   //    std::println();
   // }
   // std::println();

   // Constants indices
   auto x_idx = 2zu;
   auto y_idx = 2zu;
   auto z_idx = 2zu;

   /* ---------------------- XY 2d slice ---------------------- */
   subspan2d_t xy_span{
      &md1[0, 0, z_idx],         // Pointer to data start
      {std::extents{nx, ny},     // Slice shape
       std::array{ny * nz, 1zu}} // Slice strides
   };

   for (std::size_t i = 0; i < nx; i++) {
      for (std::size_t j = 0; j < ny; j++) {
         std::print("{:>3}, ", xy_span[i, j]);
      }
      std::println();
   }
   std::println();

   /* ---------------------- XZ 2d slice ---------------------- */
   subspan2d_t xz_span{
      &md1[0, y_idx, 0],         // Pointer to data start
      {std::extents{nx, nz},     // Slice shape
       std::array{ny * nz, 1zu}} // Slice strides
   };
   for (std::size_t i = 0; i < nx; i++) {
      for (std::size_t k = 0; k < nz; k++) {
         std::print("{:>3}, ", xz_span[i, k]);
      }
      std::println();
   }
   std::println();

   /* ---------------------- YZ 2d slice ---------------------- */
   subspan2d_t yz_span{
      &md1[x_idx, 0, 0],     // Pointer to data start
      {std::extents{ny, nz}, // Slice shape
       std::array{nz, 1zu}}  // Slice strides
   };
   for (std::size_t j = 0; j < ny; j++) {
      for (std::size_t k = 0; k < nz; k++) {
         std::print("{:>3}, ", yz_span[j, k]);
      }
      std::println();
   }
   std::println();

   /* ---------------------- x-axis 1d slice ---------------------- */
   subspan_t x_slice{
      &md1[0, y_idx, z_idx], // Pointer to data start
      {std::extents{nx},     // Slice shape
       std::array{ny * nz}}  // Slice strides
   };

   for (std::size_t i = 0; i < nx; i++) {
      std::print("{:>3}, ", x_slice[i]);
   }
   std::println();

   /* ---------------------- y-axis 1d slice ---------------------- */
   subspan_t y_slice{
      &md1[x_idx, 0, y_idx], // Pointer to data start
      {std::extents{ny},     // Slice shape
       std::array{nz}}       // Slice strides
   };

   for (std::size_t j = 0; j < ny; j++) {
      std::print("{:>3}, ", y_slice[j]);
   }
   std::println();

   /* ---------------------- z-axis 1d slice ---------------------- */
   subspan_t z_slice{
      &md1[x_idx, y_idx, 0], // Pointer to data start
      {std::extents{nz},     // Slice shape
       std::array{1zu}}      // Slice strides
   };

   for (std::size_t k = 0; k < nz; k++) {
      std::print("{:>3}, ", z_slice[k]);
   }
   std::println();
}

#endif //MDSPAN_TESTS_HPP
