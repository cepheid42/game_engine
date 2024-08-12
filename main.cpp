#include <fstream>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <ranges>
#include <sstream>
#include <functional>
#include <vector>

#define DEBUG

// #include "core/debug.h"
// #include "core/h5_wrapper.h"
// #include "core/logging.h"
// #include "core/maths.h"
// #include "core/memory.h"
// #include "core/profiling.h"
// #include "core/random.h"
// #include "core/timers.h"
// #include "core/typelist.h"
// #include "electromagnetics/em_solver.h"

#include "experimental/expression_templates.h"

void print(auto& arr) {
    for (size_t i = 0; i < arr.nx; i++) {
        for (size_t j = 0; j < arr.ny; j++) {
            std::cout << arr(i, j);
            if (j < arr.ny - 1) { std::cout << ", "; }
        }
        std::cout << '\n';
    }
    std::cout << '\n' << std::endl;
}


int main() {
    constexpr size_t nx = 5;
    constexpr size_t ny = 5;

    Matrix<int, nx, ny> A;
    Matrix<int, nx, ny> B;
    Matrix<int, nx, ny> C;

    for (size_t i = 0; i < nx; i++) {
        for (size_t j = 0; j < ny; j++) {
            A(i, j) = j + i * ny;
            B(i, j) = j;
        }
    }

    C = A + B;
    // C = diff_x(A);

    print(A);
    print(B);
    print(C);

    return 0;
}
