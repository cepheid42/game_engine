//
// Created by cepheid on 7/15/24.
//

#ifndef ARRAY_H
#define ARRAY_H

#include <cassert>

using fptype = double;

struct EMPTYARRAY {
    constexpr auto operator()(size_t...) const { return 0.0; }
};

struct Array2D {
    fptype* data;
    size_t shape[2];
    size_t size;

    Array2D() = default;

    Array2D(size_t nx, size_t ny)
    : data(new fptype[nx * ny]),
      shape{nx, ny},
      size{nx * ny}
    {}

    // ~Array() { delete[] data; }

    [[nodiscard]] auto get_index(size_t i, size_t j) const { assert(i < shape[0] and j < shape[1]); return j + (shape[1] * i); }

    // 2D indexing
    auto& operator()(const size_t i, const size_t j) { return data[get_index(i, j)]; }
    const auto& operator()(const size_t i, const size_t j) const { return data[get_index(i, j)]; }

    // 2D difference
    [[nodiscard]] auto backward_diff_x(size_t i, size_t j) const { return data[get_index(i, j)] - data[get_index(i - 1, j)]; }
    [[nodiscard]] auto backward_diff_y(size_t i, size_t j) const { return data[get_index(i, j)] - data[get_index(i, j - 1)]; }

    [[nodiscard]] auto forward_diff_x(size_t i, size_t j) const { return data[get_index(i + 1, j)] - data[get_index(i, j)]; }
    [[nodiscard]] auto forward_diff_y(size_t i, size_t j) const { return data[get_index(i, j + 1)] - data[get_index(i, j)]; }
};

struct Array3D {
    fptype* data;
    size_t shape[3];
    size_t size;

    Array3D() = default;

    Array3D(size_t nx, size_t ny, size_t nz)
    : data(new fptype[nx * ny * nz]),
      shape{nx, ny, nz},
      size{nx * ny * nz}
    {}

    // ~Array() { delete[] data; }

    [[nodiscard]] auto get_index(size_t i, size_t j, size_t k) const
    {
        assert(i < shape[0] and j < shape[1] and k < shape[2]);
        return k + (shape[2] * j) + (shape[1] * shape[2] * i);
    }

    // 3D indexing
    auto& operator()(const size_t i, const size_t j, const size_t k)
    {
        return data[get_index(i, j, k)];
    }
    const auto& operator()(const size_t i, const size_t j, const size_t k) const
    {
        return data[get_index(i, j, k)];
    }

    // 3D difference
    [[nodiscard]] auto backward_diff_x(size_t i, size_t j, size_t k) const { return data[get_index(i, j, k)] - data[get_index(i - 1, j, k)]; }
    [[nodiscard]] auto backward_diff_y(size_t i, size_t j, size_t k) const { return data[get_index(i, j, k)] - data[get_index(i, j - 1, k)]; }
    [[nodiscard]] auto backward_diff_z(size_t i, size_t j, size_t k) const { return data[get_index(i, j, k)] - data[get_index(i, j, k - 1)]; }

    [[nodiscard]] auto forward_diff_x(size_t i, size_t j, size_t k) const { return data[get_index(i + 1, j, k)] - data[get_index(i, j, k)]; }
    [[nodiscard]] auto forward_diff_y(size_t i, size_t j, size_t k) const { return data[get_index(i, j + 1, k)] - data[get_index(i, j, k)]; }
    [[nodiscard]] auto forward_diff_z(size_t i, size_t j, size_t k) const { return data[get_index(i, j, k + 1)] - data[get_index(i, j, k)]; }
};

#endif //ARRAY_H
