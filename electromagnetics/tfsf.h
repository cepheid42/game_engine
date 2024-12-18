//
// Created by cepheid on 12/16/24.
//

#ifndef TFSF_H
#define TFSF_H

#include "em_updates.h"

namespace tf::electromagnetics::sources
{

  template<typename Array>
  void print_array(const Array& arr) {
    for (size_t i = 0; i < arr.nx(); i++) {
      std::cout << arr[i] << ", ";
    }
    std::cout << std::endl;
  }

  namespace detail
  {
    inline fp_t ricker(const fp_t q) {
      constexpr auto Np = 20.0;
      constexpr auto Md = 2.5;

      const auto alpha = (M_PI * (cfl * q / Np - Md)) * (M_PI * (cfl * q / Np - Md));

      return (1.0 - 2.0 * alpha) * std::exp(-alpha);
      // return std::exp(-std::pow((q - 30.0) / 10.0, 2.0));
    }
  }

  template<typename T>
  struct TFSFSourceTM {
    using HIntegrator = FieldIntegrator1D<tf::types::Array1D<T>, FieldUpdate<Derivative::DX, Derivative::NoOp, true, size_t>>;
    using EIntegrator = FieldIntegrator1D<tf::types::Array1D<T>, FieldUpdate<Derivative::DX, Derivative::NoOp, false, size_t>>;
    using Boundary = Pml1DUpdate<EMFace::X, EMSide::Hi>;

    using pml_t = PMLData<Array1D<T>, EMFace::X, EMSide::Hi, false>;
    using empty_t = tf::types::EmptyArray1D<T>;

    static constexpr empty_t empty{};

    // todo: See about making the PML region better, or just use a dummy damping region. It could be better.
    TFSFSourceTM(size_t nx, T dt, T dx, const size_t x0_, const size_t x1_, const size_t y0_, const size_t y1_)
    : x0(x0_), x1(x1_), y0(y0_), y1(y1_),
      Einc{nx}, Ceze{nx}, Cezh{nx},
      Hinc{nx - 1}, Chyh{nx - 1}, Chye{nx - 1},
      Epml{Einc, dt, dx},
      Hpml{Hinc, dt, dx}
    {
      init_coefficients(dt, dx);
    }

    void init_coefficients(T, T);

    // This whole thing assume TMz layout with propagation in +X direction
    // todo: LSI is going to need TEy layout or something similar
    void correct_Ez(auto& Ez, const auto& Cezhy) {
      // Correct Ez @ x0
      for (size_t j = y0; j <= y1; ++j) {
        Ez(x0, j) -= Cezhy(x0 - 1, j) * Hinc[x0 - 1];
      }

      // Correct Ez @ x1
      for (size_t j = y0; j <= y1; ++j) {
        Ez(x1, j) += Cezhy(x1, j) * Hinc[x1];
      }
    }

    void correct_Hx(auto& Hx, const auto& Chxez) {
      // Correct Hx @ y0
      for (size_t i = x0; i <= x1; ++i) {
        Hx(i, y0 - 1) += Chxez(i, y0 - 1) * Einc[i];
      }

      // Correct Hx @ y1
      for (size_t i = x0; i <= x1; ++i) {
        Hx(i, y1) -= Chxez(i, y1) * Einc[i];
      }
    }


    void correct_Hy(auto& Hy, const auto& Chyez) {
      // Correct Hy @ x0
      for (size_t j = y0; j <= y1; ++j) {
        Hy(x0 - 1, j) -= Chyez(x0 - 1, j) * Einc[x0];
      }

      // Correct Hy @ x1
      for (size_t j = y0; j <= y1; ++j) {
        Hy(x1, j) += Chyez(x1, j) * Einc[x1];
      }
    }

    void apply(auto& em, const auto q) {
    // Constant X-Faces -- scattered-field nodes
    correct_Hy(em.Hy, em.Chyez);

    // Constant Y-Faces -- scattered-field nodes
    correct_Hx(em.Hx, em.Chxez);

    // Update Auxiliary 1D source
    // updateH
    HIntegrator::apply(Hinc, Einc, empty, empty, Chyh, Chye, empty, empty, {0, 0, 0, 0, 0, 0});
    Boundary::updateH(Hpml, Hinc, Einc, Chye);

    // updateE
    EIntegrator::apply(Einc, Hinc, empty, empty, Ceze, Cezh, empty, empty, {1, 1, 1, 1, 0, 0});
    Boundary::updateE(Epml, Einc, Hinc, Cezh);

    // increment Ez(0)
    Einc[0] = detail::ricker(q);

    // Constant X-Faces -- total-field nodes
    correct_Ez(em.Ez, em.Cezhy);

    // Constant Y-Faces -- total-field nodes
    // Nothing Here
    }

    size_t x0, x1, y0, y1;

    Array1D<T> Einc;
    Array1D<T> Ceze;
    Array1D<T> Cezh;

    Array1D<T> Hinc;
    Array1D<T> Chyh;
    Array1D<T> Chye;

    pml_t Epml;
    pml_t Hpml;
  };

  template <typename T>
  void TFSFSourceTM<T>::init_coefficients(T dt, T dx) {
    // constexpr auto eps0 = 8.854187812813e-12;
    // constexpr auto mu0 = 1.2566370621219e-6;
    // constexpr auto sigma = 0.0;
    //
    // // half dt for H field, since it's split into two steps
    // const auto hc = dt / (mu0 * dx);
    //
    // const auto e_num = dt / (eps0 * dx);
    // const auto alpha = (sigma * dt) / (2.0 * eps0);
    // const auto ec = (1.0 - alpha) / (1.0 + alpha);
    // const auto eh = e_num / (1.0 + alpha);

    constexpr auto eta0 = 377.0;

    init_coeff(Ceze, 1.0);
    init_coeff(Cezh, cfl * eta0);

    init_coeff(Chyh, 1.0);
    init_coeff(Chye, cfl / eta0);
  }
} // end namespace tf::electromagnetics::sources
#endif //TFSF_H
