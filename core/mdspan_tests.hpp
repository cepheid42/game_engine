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

void test_periodic() {
      constexpr auto nx = 10zu;
   constexpr auto ny = 10zu;
   constexpr auto nz = 10zu;
   constexpr auto nhalo = 2zu;

   using mdspan_t = std::mdspan<
      std::size_t,
      std::dextents<std::size_t, 3>
   >;

   using strided_slice_t = std::mdspan<std::size_t, std::dextents<std::size_t, 3>, std::layout_stride>;

   std::vector<std::size_t> data(nx * ny * nz);
   mdspan_t test{data.data(), std::extents{nx, ny, nz}};

   // for (auto i = 0zu; i < nhalo; i++) {
   //    for (auto j = 0zu; j < test.extent(1); j++) {
   //       for (auto k = 0zu; k < test.extent(2); k++) {
   //          test[i, j, k] = i + 1;
   //       }
   //    }
   // }
   //
   // for (auto i = nx - nhalo; i < test.extent(0); i++) {
   //    for (auto j = 0zu; j < test.extent(1); j++) {
   //       for (auto k = 0zu; k < test.extent(2); k++) {
   //          test[i, j, k] = i + 1;
   //       }
   //    }
   // }

   // for (auto i = 0zu; i < test.extent(0); i++) {
   //    for (auto j = 0zu; j < nhalo; j++) {
   //       for (auto k = 0zu; k < test.extent(2); k++) {
   //          test[i, j, k] = j + 1;
   //       }
   //    }
   // }
   //
   // for (auto i = 0zu; i < test.extent(0); i++) {
   //    for (auto j = ny - nhalo; j < test.extent(1); j++) {
   //       for (auto k = 0zu; k < test.extent(2); k++) {
   //          test[i, j, k] = j + 1;
   //       }
   //    }
   // }

   for (auto i = 0zu; i < test.extent(0); i++) {
      for (auto j = 0zu; j < test.extent(1); j++) {
         for (auto k = 0zu; k < nhalo; k++) {
            test[i, j, k] = k + 1;
         }
      }
   }

   for (auto i = 0zu; i < test.extent(0); i++) {
      for (auto j = 0zu; j < test.extent(1); j++) {
         for (auto k = nz - nhalo; k < test.extent(2); k++) {
            test[i, j, k] = k + 1;
         }
      }
   }

   // stepfill_mdarray(test, 0);
   // print_mdarray(test);

   // // X dir
   // strided_slice_t  first{&test[0, 0, 0], {std::extents{1, ny, nz}, std::array{(nx - 2 * nhalo) * ny * nz, nz, 1zu}}}; // first plane
   // strided_slice_t second{&test[1, 0, 0], {std::extents{1, ny, nz}, std::array{(nx - 2 * nhalo) * ny * nz, nz, 1zu}}}; // second plane
   // strided_slice_t  third{&test[2, 0, 0], {std::extents{1, ny, nz}, std::array{(nx - 2 * nhalo) * ny * nz, nz, 1zu}}}; // second last plane
   // strided_slice_t fourth{&test[3, 0, 0], {std::extents{1, ny, nz}, std::array{(nx - 2 * nhalo) * ny * nz, nz, 1zu}}}; // last plane

   // // Y dir
   // strided_slice_t  first{&test[0, 0, 0], {std::extents{nx, 1, nz}, std::array{ny * nz, (ny - 2 * nhalo) * nz, 1zu}}};
   // strided_slice_t second{&test[0, 1, 0], {std::extents{nx, 2, nz}, std::array{ny * nz, (ny - 2 * nhalo) * nz, 1zu}}};
   // strided_slice_t  third{&test[0, 2, 0], {std::extents{nx, 2, nz}, std::array{ny * nz, (ny - 2 * nhalo) * nz, 1zu}}};
   // strided_slice_t fourth{&test[0, 3, 0], {std::extents{nx, 2, nz}, std::array{ny * nz, (ny - 2 * nhalo) * nz, 1zu}}};

   // Z dir
   strided_slice_t  first{&test[0, 0, 0], {std::extents{nx, ny, 1}, std::array{ny * nz, nz, nz - 2 * nhalo}}};
   strided_slice_t second{&test[0, 0, 1], {std::extents{nx, ny, 1}, std::array{ny * nz, nz, nz - 2 * nhalo}}};
   strided_slice_t  third{&test[0, 0, 2], {std::extents{nx, ny, 1}, std::array{ny * nz, nz, nz - 2 * nhalo}}};
   strided_slice_t fourth{&test[0, 0, 3], {std::extents{nx, ny, 1}, std::array{ny * nz, nz, nz - 2 * nhalo}}};

   for (auto i = 0u; i < first.extent(0); i++) {
      for (auto j = 0u; j < first.extent(1); j++) {
         for (auto k = 0u; k < first.extent(2); k++) {
            // first[i + 1, j, k] = first[i, j, k];
            // second[i + 1, j, k] = second[i, j, k];
            // third[i, j, k] = third[i + 1, j, k];
            // fourth[i, j, k] = fourth[i + 1, j, k];

            // first[i, j + 1, k] = first[i, j, k];
            // second[i, j + 1, k] = second[i, j, k];
            // third[i, j, k] = third[i, j + 1, k];
            // fourth[i, j, k] = fourth[i, j + 1, k];

            first[i, j, k + 1] = first[i, j, k];
            second[i, j, k + 1] = second[i, j, k];
            third[i, j, k] = third[i, j, k + 1];
            fourth[i, j, k] = fourth[i, j, k + 1];
         }
      }
   }

   print_mdarray(test);
}
#endif //MDSPAN_TESTS_HPP
