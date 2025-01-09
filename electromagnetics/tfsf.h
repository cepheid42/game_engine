//
// Created by cepheid on 12/16/24.
//

#ifndef TFSF_H
#define TFSF_H

// todo: TFSF was only implemented for TMz mode, no 3D version yet
//       so this file is useless right now.

//#include <bc_data.h>
//#include <boundaries.h>
//
//#include "em_updates.h"
//
//namespace tf::electromagnetics::sources
//{
//
//  template<typename Array>
//  void print_array(const Array& arr) {
//    for (size_t i = 0; i < arr.nx(); i++) {
//      std::cout << arr[i] << ", ";
//    }
//    std::cout << std::endl;
//  }
//
//  namespace detail
//  {
//    inline fp_t ricker(const fp_t q) {
//      constexpr auto Np = 20.0;
//      constexpr auto Md = 2.0;
//
//      const auto alpha = (M_PI * (cfl * q / Np - Md)) * (M_PI * (cfl * q / Np - Md));
//
//      return (1.0 - 2.0 * alpha) * std::exp(-alpha);
//      // return std::exp(-std::pow((q - 30.0) / 10.0, 2.0));
//    }
//  }
//
//  template<typename T>
//  struct TFSFSourceTM {
//    using HIntegrator = FieldIntegrator3D<tf::types::Array3D<T>, FieldUpdate<Derivative::DX, Derivative::NoOp, true, size_t>>;
//    using EIntegrator = FieldIntegrator3D<tf::types::Array3D<T>, FieldUpdate<Derivative::DX, Derivative::NoOp, false, size_t>>;
//    using empty_t = tf::types::EmptyArray3D<T>;
//
//    static constexpr empty_t empty{};
//
//    TFSFSourceTM(const size_t nx, const T dt, const T dx, const size_t x0_, const size_t x1_, const size_t y0_, const size_t y1_)
//    : x0(x0_), x1(x1_), y0(y0_), y1(y1_),
//      Einc{nx + 20}, Ceze{nx + 20}, Cezh{nx + 20},
//      Hinc{nx - 1 + 20}, Chyh{nx - 1 + 20}, Chye{nx - 1 + 20}
//    {
//      // todo: Adds an additional 20 nodes. Need to make sure this works when nPML == 0
//      //       since below the NLOSS region is set to be 30 nodes
//      init_coefficients(dt, dx);
//    }
//
//    void init_coefficients(T, T);
//
//    // todo: This whole thing assume TMz layout with propagation in +X direction
//    //       LSI is going to need TEy layout or something similar
//    void correct_Ez(auto& Ez, const auto& Cezhy) {
//      // Correct Ez @ x0
//      for (size_t j = y0; j <= y1; ++j) {
//        // todo: make sure offsetting by the PML depth is a good idea later on
//        Ez(x0, j) -= Cezhy(x0, j) * Hinc[x0 - 1 - nPml];
//      }
//
//      // Correct Ez @ x1
//      for (size_t j = y0; j <= y1; ++j) {
//        Ez(x1, j) += Cezhy(x1, j) * Hinc[x1 - nPml];
//      }
//    }
//
//    void correct_Hx(auto& Hx, const auto& Chxez) {
//      // Correct Hx @ y0
//      for (size_t i = x0; i <= x1; ++i) {
//        Hx(i, y0 - 1, 0) += 2.0 * Chxez(i, y0 - 1) * Einc[i - nPml];
//      }
//
//      // Correct Hx @ y1
//      for (size_t i = x0; i <= x1; ++i) {
//        Hx(i, y1) -= 2.0 * Chxez(i, y1) * Einc[i - nPml];
//      }
//    }
//
//    void correct_Hy(auto& Hy, const auto& Chyez) {
//      // Correct Hy @ x0
//      for (size_t j = y0; j <= y1; ++j) {
//        Hy(x0 - 1, j) -= 2.0 * Chyez(x0 - 1, j) * Einc[x0 - nPml];
//      }
//
//      // Correct Hy @ x1
//      for (size_t j = y0; j <= y1; ++j) {
//        Hy(x1, j) += 2.0 * Chyez(x1, j) * Einc[x1 - nPml];
//      }
//    }
//
//    void apply(auto& em, const auto q) {
//    // Constant X-Faces -- scattered-field nodes
//    correct_Hy(em.Hy, em.Chyez);
//
//    // Constant Y-Faces -- scattered-field nodes
//    correct_Hx(em.Hx, em.Chxez);
//
//    // Update Auxiliary 1D source
//    // updateH
//    HIntegrator::apply(Hinc, Einc, empty, empty, Chyh, Chye, empty, empty, {0, 0, 0, 0, 0, 0});
//    // updateE
//    EIntegrator::apply(Einc, Hinc, empty, empty, Ceze, Cezh, empty, empty, {1, 1, 1, 1, 0, 0});
//
//    // increment Ez(0)
//    // todo: Figure out how to pass different functions in as template params
//    Einc[0] = detail::ricker(q);
//
//    // Constant X-Faces -- total-field nodes
//    correct_Ez(em.Ez, em.Cezhy);
//
//    // Constant Y-Faces -- total-field nodes
//    // Nothing Here
//    }
//
//    size_t x0, x1, y0, y1;
//
//    tf::types::Array3D<T> Einc;
//    tf::types::Array3D<T> Ceze;
//    tf::types::Array3D<T> Cezh;
//
//    tf::types::Array3D<T> Hinc;
//    tf::types::Array3D<T> Chyh;
//    tf::types::Array3D<T> Chye;
//  };
//
//  template <typename T>
//  void TFSFSourceTM<T>::init_coefficients(T dt, T dx) {
//    // todo: these constants are in a different header
//    constexpr auto eps0 = 8.854187812813e-12;
//    constexpr auto mu0 = 1.2566370621219e-6;
//
//    // todo: Need to incorporate the grids eps_r and mu_r values into these coeffs
//    const auto hc = dt / (mu0 * dx);
//    const auto eh = dt / (eps0 * dx);
//
//    //
//    constexpr size_t NLOSS = 20;
//    constexpr double MAX_LOSS = 0.35;
//
//    const auto sizex = Einc.nx();
//
//    for (size_t i = 0; i < sizex - 1; ++i) {
//      if (i < sizex - 1 - NLOSS) {
//        Ceze[i] = 1.0;
//        Cezh[i] = eh;
//        Chyh[i] = 1.0;
//        Chye[i] = hc;
//      } else {
//        auto depth = static_cast<double>(i - (sizex - 1 - NLOSS)) + 0.5;
//        auto loss = MAX_LOSS * std::pow(depth / NLOSS, 3.0);
//        Ceze[i] = (1.0 - loss) / (1.0 + loss);
//        Cezh[i] = eh / (1.0 + loss);
//
//        depth += 0.5;
//        loss = MAX_LOSS * std::pow(depth / NLOSS, 3.0);
//        Chyh[i] = (1.0 - loss) / (1.0 + loss);
//        Chye[i] = hc / (1.0 + loss);
//      }
//    }
//  }
//
//} // end namespace tf::electromagnetics::sources
#endif //TFSF_H
