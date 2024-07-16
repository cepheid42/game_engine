#include <fstream>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <functional>
#include <vector>
#include <algorithm>
#include <numeric>

#define DEBUG

#include "core/dbg.h"
#include "core/typelist.h"
#include "core/timer.h"
#include "electromagnetics/array.h"
#include "electromagnetics/em_solver.h"
#include "electromagnetics/regular_em.h"

using fptype = double;

inline constexpr double pi{3.141592653589793};
inline constexpr double eta0{376.73031366686992}; // Ohms
// inline constexpr double cfl{1.0 / 1.414213562373095}; // 2D = 1 / sqrt(2)
inline constexpr double cfl{1.0 / 1.732050807568877};    // 3D = 1 / sqrt(3)


auto gauss_deriv(const size_t n) {
    const auto q = static_cast<double>(n);
    const auto c = -(2.0 * (q - 40.0) / (10.0 * 10.0));
    const auto a = (q - 40.0) / 10.0;
    const auto aa = a * a;
    const auto result = c * std::exp(-aa);
    return result;
}

struct FieldData2D {
    FieldData2D(size_t nx, size_t ny)
    : Ez(nx, ny),
      Hx(nx, ny - 1),
      Hy(nx - 1, ny)
    {}
    
    Array2D Ez;
    Array2D Hx;
    Array2D Hy;
    EMPTYARRAY Ex{};
    EMPTYARRAY Jx{};
    EMPTYARRAY Jz{};
    EMPTYARRAY Ey{};
    EMPTYARRAY Hz{};
    EMPTYARRAY Jy{};
};

struct CoeffData2D
{
    CoeffData2D(size_t nx, size_t ny)
    : ceze(nx, ny),
      cezh(nx, ny),
      chxh(nx, ny - 1),
      chxe(nx, ny - 1),
      chyh(nx - 1, ny),
      chye(nx - 1, ny)
    {}
    Array2D ceze;
    Array2D cezh;
    Array2D chxh;
    Array2D chxe;
    Array2D chyh;
    Array2D chye;
    EMPTYARRAY cexe{};
    EMPTYARRAY cexh{};
    EMPTYARRAY cjx{};
    EMPTYARRAY cjz{};
    EMPTYARRAY ceye{};
    EMPTYARRAY ceyh{};
    EMPTYARRAY chzh{};
    EMPTYARRAY chze{};
    EMPTYARRAY cjy{};
};
void set_coefficients2D(CoeffData2D& c)
{
    for (size_t i = 0; i < c.ceze.shape[0]; i++) {
        for (size_t j = 0; j < c.ceze.shape[1]; j++) {
            c.ceze(i, j) = 1.0;
            c.cezh(i, j) = cfl * eta0;
        }
    }

    for (size_t i = 0; i < c.chxh.shape[0]; i++) {
        for (size_t j = 0; j < c.chxh.shape[1]; j++) {
            c.chxh(i, j) = 1.0;
            c.chxe(i, j) = cfl / eta0;
        }
    }

    for (size_t i = 0; i < c.chyh.shape[0]; i++) {
        for (size_t j = 0; j < c.chyh.shape[1]; j++) {
            c.chyh(i, j) = 1.0;
            c.chye(i, j) = cfl / eta0;
        }
    }
}
void save_array(const Array2D& a, const std::string& path, const size_t step)
{
    std::string count_padded = std::to_string(step);
    count_padded.insert(count_padded.begin(), 8 - count_padded.length(), '0');

    const auto filename = path + '_' + count_padded + ".csv";

    std::ofstream output_file(filename, std::ios::trunc);
    for (size_t i = 0; i < a.shape[0]; i++) {
        for (size_t j = 0; j < a.shape[1]; j++) {
            output_file << a(i, j);
            if (j < a.shape[1] - 1) {
                output_file << ", ";
            }
        }
        output_file << '\n';
    }
    output_file.close();
}

template<typename Array>
auto calculate_rms2D(const Array& fancy, const Array& regular, const std::string& name) {
    assert(fancy.shape[0] == regular.shape[0]);
    assert(fancy.shape[1] == regular.shape[1]);
    std::vector<fptype> result(fancy.size);
    for (size_t i = 0; i < fancy.shape[0]; i++) {
        for (size_t j = 0; j < fancy.shape[1]; j++) {
            const auto r = fancy(i, j) - regular(i, j);
            result.push_back(r * r);
        }
    }
    auto lse = std::accumulate(result.begin(), result.end(), 0.0);
    lse = std::sqrt(lse / static_cast<fptype>(fancy.size));
    std::cout << name << ": " << std::scientific << lse << std::endl;
}

// struct FieldData3D {
//     FieldData3D(size_t nx, size_t ny, size_t nz)
//     : Ex(nx - 1, ny, nz),
//       Ey(nx, ny - 1, nz),
//       Ez(nx, ny, nz - 1),
//       Hx(nx, ny - 1, nz - 1),
//       Hy(nx - 1, ny, nz - 1),
//       Hz(nx - 1, ny - 1, nz)
//     {}
//     Array3D Ex;
//     Array3D Ey;
//     Array3D Ez;
//     Array3D Hx;
//     Array3D Hy;
//     Array3D Hz;
//     EMPTYARRAY Jx{};
//     EMPTYARRAY Jz{};
//     EMPTYARRAY Jy{};
// };
//
//
// struct CoeffData3D {
//     CoeffData3D(size_t nx, size_t ny, size_t nz)
//     : cexe(nx - 1, ny, nz),
//       cexh(nx - 1, ny, nz),
//       ceye(nx, ny - 1, nz),
//       ceyh(nx, ny - 1, nz),
//       ceze(nx, ny, nz - 1),
//       cezh(nx, ny, nz - 1),
//       chxh(nx, ny - 1, nz - 1),
//       chxe(nx, ny - 1, nz - 1),
//       chyh(nx - 1, ny, nz - 1),
//       chye(nx - 1, ny, nz - 1),
//       chzh(nx - 1, ny - 1, nz),
//       chze(nx - 1, ny - 1, nz)
//     {}
//     Array3D cexe;
//     Array3D cexh;
//     Array3D ceye;
//     Array3D ceyh;
//     Array3D ceze;
//     Array3D cezh;
//     Array3D chxh;
//     Array3D chxe;
//     Array3D chyh;
//     Array3D chye;
//     Array3D chzh;
//     Array3D chze;
//     EMPTYARRAY cjx{};
//     EMPTYARRAY cjz{};
//     EMPTYARRAY cjy{};
// };
// void set_coefficients3D(CoeffData3D& c) {
//     for (size_t i = 0; i < c.cexe.shape[0]; i++) {
//         for (size_t j = 0; j < c.cexe.shape[1]; j++) {
//             for (size_t k = 0; k < c.cexe.shape[2]; k++) {
//                 c.cexe(i, j, k) = 1.0;
//                 c.cexh(i, j, k) = cfl * eta0;
//             }
//         }
//     }
//     for (size_t i = 0; i < c.ceye.shape[0]; i++) {
//         for (size_t j = 0; j < c.ceye.shape[1]; j++) {
//             for (size_t k = 0; k < c.ceye.shape[2]; k++) {
//                 c.ceye(i, j, k) = 1.0;
//                 c.ceyh(i, j, k) = cfl * eta0;
//             }
//         }
//     }
//     for (size_t i = 0; i < c.ceze.shape[0]; i++) {
//         for (size_t j = 0; j < c.ceze.shape[1]; j++) {
//             for (size_t k = 0; k < c.ceze.shape[2]; k++) {
//                 c.ceze(i, j, k) = 1.0;
//                 c.cezh(i, j, k) = cfl * eta0;
//             }
//         }
//     }
//     for (size_t i = 0; i < c.chxh.shape[0]; i++) {
//         for (size_t j = 0; j < c.chxh.shape[1]; j++) {
//             for (size_t k = 0; k < c.chxh.shape[2]; k++) {
//                 c.chxh(i, j, k) = 1.0;
//                 c.chxe(i, j, k) = cfl / eta0;
//             }
//         }
//     }
//     for (size_t i = 0; i < c.chyh.shape[0]; i++) {
//         for (size_t j = 0; j < c.chyh.shape[1]; j++) {
//             for (size_t k = 0; k < c.chyh.shape[2]; k++) {
//                 c.chyh(i, j, k) = 1.0;
//                 c.chye(i, j, k) = cfl / eta0;
//             }
//         }
//     }
//     for (size_t i = 0; i < c.chzh.shape[0]; i++) {
//         for (size_t j = 0; j < c.chzh.shape[1]; j++) {
//             for (size_t k = 0; k < c.chzh.shape[2]; k++) {
//                 c.chzh(i, j, k) = 1.0;
//                 c.chze(i, j, k) = cfl / eta0;
//             }
//         }
//     }
// }
//
// void save_array(const Array3D& a, const std::string& path, const size_t step) {
//     std::string count_padded = std::to_string(step);
//     count_padded.insert(count_padded.begin(), 8 - count_padded.length(), '0');
//     const auto filename = path + '_' + count_padded + ".csv";
//     std::ofstream output_file(filename, std::ios::trunc);
//     for (size_t i = 0; i < a.shape[0]; i++) {
//         for (size_t j = 0; j < a.shape[1]; j++) {
//             for (size_t k = 0; k < a.shape[2]; k++) {
//                 output_file << a(i, j, k);
//                 if (k < a.shape[2] - 1) {
//                     output_file << ", ";
//                 }
//             }
//             output_file << '\n';
//         }
//     }
//     output_file.close();
// }
//
// template<typename Array>
// auto calculate_rms3D(const Array& fancy, const Array& regular, const std::string& name) {
//     assert(fancy.shape[0] == regular.shape[0]);
//     assert(fancy.shape[1] == regular.shape[1]);
//     assert(fancy.shape[2] == regular.shape[2]);
//     std::vector<fptype> result(fancy.size);
//     for (size_t i = 0; i < fancy.shape[0]; i++) {
//         for (size_t j = 0; j < fancy.shape[1]; j++) {
//             for (size_t k = 0; k < fancy.shape[2]; k++) {
//                 const auto r = fancy(i, j, k) - regular(i, j, k);
//                 result.push_back(r * r);
//             }
//         }
//     }
//     auto lse = std::accumulate(result.begin(), result.end(), 0.0);
//     lse = std::sqrt(lse / static_cast<fptype>(fancy.size));
//     std::cout << name << ": " << std::scientific << lse << std::endl;
// }

int main()
{
    constexpr size_t nx = 100;
    constexpr size_t ny = 100;
    constexpr size_t nz = 100;

    constexpr size_t nt = 200;
    
    // constexpr auto save_step = 2lu;
    // size_t file_count = 0;

    using ARRAY = Array2D;
    using EMPTY = EMPTYARRAY;
    
    //                            E      B1     B2     J      C1     C2     C3
    using EzUpdateTL = TypeList<ARRAY, ARRAY, ARRAY, EMPTY, ARRAY, ARRAY, EMPTY>;
    //                            H      E1     E2     M      C1     C2     C3
    using HxUpdateTL = TypeList<ARRAY, EMPTY, ARRAY, EMPTY, ARRAY, ARRAY, EMPTY>;
    using HyUpdateTL = TypeList<ARRAY, ARRAY, EMPTY, EMPTY, ARRAY, ARRAY, EMPTY>;
    using ExIntegrator = FieldIntegratorNull;
    using EyIntegrator = FieldIntegratorNull;
    using EzIntegrator = FieldIntegrator2D<EzUpdateTL, FieldType::E, Derivative::DX, Derivative::DY>;
    using HxIntegrator = FieldIntegrator2D<HxUpdateTL, FieldType::H, Derivative::NoOp, Derivative::DY>;
    using HyIntegrator = FieldIntegrator2D<HyUpdateTL, FieldType::H, Derivative::DX, Derivative::NoOp>;
    using HzIntegrator = FieldIntegratorNull;
    using EMSolver = Electromagnetics<ExIntegrator, EyIntegrator, EzIntegrator, HxIntegrator, HyIntegrator, HzIntegrator>;
    FieldData2D ff(nx, ny);
    FieldData2D rf(nx, ny);
    CoeffData2D coeffs(nx, ny);
    set_coefficients2D(coeffs);
    for (size_t n = 0; n < nt; n++) {
        std::cout << "Step " << n << std::endl;
        const auto src = gauss_deriv(n);
        EMSolver::updateB(ff, coeffs);
        EMSolver::updateE(ff, coeffs);
        ff.Ez(50, 50) += src;
        update_Hx2D(rf.Hx, rf.Ez, coeffs.chxh, coeffs.chxe);
        update_Hy2D(rf.Hy, rf.Ez, coeffs.chyh, coeffs.chye);
        update_Ez2D(rf.Ez, rf.Hy, rf.Hx, coeffs.ceze, coeffs.cezh);
        rf.Ez(50, 50) += src;
        // if (n % save_step == 0) {
        //     save_array(fields.Ez, "/home/cepheid/TriForce/game_engine/data/ez", file_count);
        //     file_count++;
        // }
    }
    calculate_rms2D(ff.Ez, rf.Ez, "Ez");
    calculate_rms2D(ff.Hx, rf.Hx, "Hx");
    calculate_rms2D(ff.Hy, rf.Hy, "Hy");

    // //                            E      H1    H2      J      C1     C2     C3
    // using ExUpdateTL = TypeList<ARRAY, ARRAY, ARRAY, EMPTY, ARRAY, ARRAY, EMPTY>;
    // using EyUpdateTL = TypeList<ARRAY, ARRAY, ARRAY, EMPTY, ARRAY, ARRAY, EMPTY>;
    // using EzUpdateTL = TypeList<ARRAY, ARRAY, ARRAY, EMPTY, ARRAY, ARRAY, EMPTY>;
    // //                            H      E1    E2      M      C1     C2     C3
    // using HxUpdateTL = TypeList<ARRAY, ARRAY, ARRAY, EMPTY, ARRAY, ARRAY, EMPTY>;
    // using HyUpdateTL = TypeList<ARRAY, ARRAY, ARRAY, EMPTY, ARRAY, ARRAY, EMPTY>;
    // using HzUpdateTL = TypeList<ARRAY, ARRAY, ARRAY, EMPTY, ARRAY, ARRAY, EMPTY>;
    //
    // using ExIntegrator = FieldIntegrator3D<ExUpdateTL, FieldType::E, Derivative::DY, Derivative::DZ>;
    // using EyIntegrator = FieldIntegrator3D<EyUpdateTL, FieldType::E, Derivative::DZ, Derivative::DX>;
    // using EzIntegrator = FieldIntegrator3D<EzUpdateTL, FieldType::E, Derivative::DX, Derivative::DY>;
    //
    // using HxIntegrator = FieldIntegrator3D<HxUpdateTL, FieldType::H, Derivative::DZ, Derivative::DY>;
    // using HyIntegrator = FieldIntegrator3D<HyUpdateTL, FieldType::H, Derivative::DX, Derivative::DZ>;
    // using HzIntegrator = FieldIntegrator3D<HzUpdateTL, FieldType::H, Derivative::DY, Derivative::DX>;
    //
    // using EMSolver = Electromagnetics<ExIntegrator, EyIntegrator, EzIntegrator, HxIntegrator, HyIntegrator, HzIntegrator>;
    // FieldData3D ff(nx, ny, nz);
    // FieldData3D rf(nx, ny, nz);
    // CoeffData3D coeffs(nx, ny, nz);
    // set_coefficients3D(coeffs);
    // for (size_t n = 0; n < nt; n++) {
    //     std::cout << "Step " << n << std::endl;
    //     const auto src = gauss_deriv(n);
    //     EMSolver::updateB(ff, coeffs);
    //     EMSolver::updateE(ff, coeffs);
    //     ff.Ez(50, 50, 50) += src;
    //     update_Hx3D(rf.Hx, rf.Ey, rf.Ez, coeffs.chxh, coeffs.chxe);
    //     update_Hy3D(rf.Hy, rf.Ez, rf.Ex, coeffs.chyh, coeffs.chye);
    //     update_Hz3D(rf.Hz, rf.Ex, rf.Ey, coeffs.chzh, coeffs.chze);
    //     update_Ex3D(rf.Ex, rf.Hz, rf.Hy, coeffs.cexe, coeffs.cexh);
    //     update_Ey3D(rf.Ey, rf.Hx, rf.Hz, coeffs.ceye, coeffs.ceyh);
    //     update_Ez3D(rf.Ez, rf.Hy, rf.Hx, coeffs.ceze, coeffs.cezh);
    //     rf.Ez(50, 50, 50) += src;
    //     // if (n % save_step == 0) {
    //     //     save_array(fields.Ez, "/home/cepheid/TriForce/game_engine/data/ez", file_count);
    //     //     file_count++;
    //     // }
    // }
    // calculate_rms3D(ff.Ex, rf.Ex, "Ex");
    // calculate_rms3D(ff.Ey, rf.Ey, "Ey");
    // calculate_rms3D(ff.Ez, rf.Ez, "Ez");
    // calculate_rms3D(ff.Hx, rf.Hx, "Hx");
    // calculate_rms3D(ff.Hy, rf.Hy, "Hy");
    // calculate_rms3D(ff.Hz, rf.Hz, "Hz");
    
    return 0;
}
