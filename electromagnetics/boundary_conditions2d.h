//
// Created by cepheid on 11/21/23.
//

#ifndef TRIFORCE_EM_BOUNDARY_CONDITIONS_H
#define TRIFORCE_EM_BOUNDARY_CONDITIONS_H

#include <electromagnetics_traits.h>

#include <algorithm>
#include <vector>

#include "array2d.h"
#include "constants.h"
#include "electromagnetics.param"
#include "electromagnetics_traits.h"
#include "mesh2d.h"
#include "program.param"
#include "typelist.h"
#include "utilities.h"
#include "emdata/emdata.h"

namespace electromagnetics::boundary_conditions
{
  enum class FieldSide : size_t { LO, HI };

  /******************************************/
  /********** PML Helper Functions **********/
  // [Taflove 7.66]
  constexpr fptype sigma_optimal(fptype dx) {
    // todo: needs relative permittivity/permeability
    //       and incorporate dz
    return (0.8 * (kPMLGradingOrder + 1.0)) / (constants::eta0 * dx);
  }

  // [Taflove 7.60a,b]
  inline std::vector<fptype> calculate_sigma(const std::vector<fptype>& d, fptype delta) {
    std::vector<fptype> sigma(d);
    auto sigma_max = sigma_optimal(delta);

    for (auto& x: sigma) {
      x = sigma_max * std::pow(x, kPMLGradingOrder);
    }

    return sigma;
  }

  // [Taflove 7.60a,b]
  inline std::vector<fptype> calculate_kappa(const std::vector<fptype>& d, fptype kappa_max) {
    std::vector<fptype> kappa(d);
    for (auto& x: kappa) {
      x = 1.0 + (kappa_max - 1.0) * std::pow(x, kPMLGradingOrder);
    }
    return kappa;
  }

  // [Taflove 7.79]
  inline std::vector<fptype> calculate_alpha(const std::vector<fptype>& d, fptype alpha_max) {
    std::vector<fptype> alpha(d);
    for (auto& x: alpha) {
      x = alpha_max * std::pow(x, kPMLAlphaOrder);
    }
    return alpha;
  }

  // [Taflove 7.102 and 7.99]
  inline void calculate_pml_coefs(std::vector<fptype>& b,
                                  std::vector<fptype>& c,
                                  const std::vector<fptype>& sigma,
                                  const std::vector<fptype>& kappa,
                                  const std::vector<fptype>& alpha,
                                  fptype dt,
                                  fptype dx,
                                  size_t d_pml)
  {
    const auto c1 = -dt / constants::epsilon0;
    for (size_t i = 0; i < d_pml; ++i) {
      b[i] = std::exp(c1 * ((sigma[i] / kappa[i]) + alpha[i]));
      c[i] = (sigma[i] * (b[i] - 1.0)) / (dx * kappa[i] * (sigma[i] + (kappa[i] * alpha[i])));
    }
  }

  /************************************************/
  /********** Boundary Condition Classes **********/
  template<BoundaryType B, FieldAxes A>
  struct Boundary {
    static constexpr auto boundary_v = B;
    static constexpr auto axes_v = A;
  };

  template<FieldAxes A>
  struct Boundary<BoundaryType::Periodic, A> {
    static constexpr auto boundary_v = BoundaryType::Periodic;
    static constexpr auto axes_v = A;
  };

  template<FieldAxes A>
  struct Boundary<BoundaryType::PML, A> {
    static constexpr auto boundary_v = BoundaryType::PML;
    static constexpr auto axes_v = A;

    explicit Boundary(size_t nx, size_t nz);

    Array2D<fptype> psiE;
    Array2D<fptype> psiB;
    std::vector<fptype> bE;
    std::vector<fptype> cE;
    std::vector<fptype> bB;
    std::vector<fptype> cB;
  };

  template<>
  inline Boundary<BoundaryType::PML, FieldAxes::X>::Boundary(size_t nx, size_t nz)
  : psiE(nx, nz + 1),
    psiB(nx, nz + 1),
    bE(nx),
    cE(nx),
    bB(nx),
    cB(nx)
  {}

  template<>
  inline Boundary<BoundaryType::PML, FieldAxes::Z>::Boundary(size_t nx, size_t nz)
  : psiE(nx + 1, nz),
    psiB(nx + 1, nz),
    bE(nz),
    cE(nz),
    bB(nz),
    cB(nz)
  {}

  /**********************************************/
  /********** Boundary Factory Methods **********/
  template<BoundaryType T, FieldAxes A, FieldSide S>
  Boundary<T, A> createBoundary(const Mesh2D& mesh, size_t nx, size_t nz, fptype dt, size_t d_pml);

  template<>
  Boundary<BoundaryType::None, FieldAxes::X>
  inline createBoundary<BoundaryType::None, FieldAxes::X, FieldSide::LO>(
    [[maybe_unused]]const Mesh2D& mesh,
    [[maybe_unused]] size_t nx,
    [[maybe_unused]] size_t nz,
    [[maybe_unused]] fptype dt,
    [[maybe_unused]] size_t d_pml)
  {
    return Boundary<BoundaryType::None, FieldAxes::X>{};
  }

  template<>
  Boundary<BoundaryType::None, FieldAxes::X>
  inline createBoundary<BoundaryType::None, FieldAxes::X, FieldSide::HI>(
    [[maybe_unused]]const Mesh2D& mesh,
    [[maybe_unused]] size_t nx,
    [[maybe_unused]] size_t nz,
    [[maybe_unused]] fptype dt,
    [[maybe_unused]] size_t d_pml)
  {
    return Boundary<BoundaryType::None, FieldAxes::X>{};
  }

  template<>
  Boundary<BoundaryType::None, FieldAxes::Z>
  inline createBoundary<BoundaryType::None, FieldAxes::Z, FieldSide::LO>(
    [[maybe_unused]]const Mesh2D& mesh,
    [[maybe_unused]] size_t nx,
    [[maybe_unused]] size_t nz,
    [[maybe_unused]] fptype dt,
    [[maybe_unused]] size_t d_pml)
  {
    return Boundary<BoundaryType::None, FieldAxes::Z>{};
  }

  template<>
  Boundary<BoundaryType::None, FieldAxes::Z>
  inline createBoundary<BoundaryType::None, FieldAxes::Z, FieldSide::HI>(
    [[maybe_unused]]const Mesh2D& mesh,
    [[maybe_unused]] size_t nx,
    [[maybe_unused]] size_t nz,
    [[maybe_unused]] fptype dt,
    [[maybe_unused]] size_t d_pml)
  {
    return Boundary<BoundaryType::None, FieldAxes::Z>{};
  }

  template<>
  Boundary<BoundaryType::Periodic, FieldAxes::Z>
  inline createBoundary<BoundaryType::Periodic, FieldAxes::Z, FieldSide::LO>(
    [[maybe_unused]]const Mesh2D& mesh,
    [[maybe_unused]] size_t nx,
    [[maybe_unused]] size_t nz,
    [[maybe_unused]] fptype dt,
    [[maybe_unused]] size_t d_pml)
  {
    return Boundary<BoundaryType::Periodic, FieldAxes::Z>{};
  }

  template<>
  Boundary<BoundaryType::Periodic, FieldAxes::Z>
  inline createBoundary<BoundaryType::Periodic, FieldAxes::Z, FieldSide::HI>(
    [[maybe_unused]]const Mesh2D& mesh,
    [[maybe_unused]] size_t nx,
    [[maybe_unused]] size_t nz,
    [[maybe_unused]] fptype dt,
    [[maybe_unused]] size_t d_pml)
  {
    return Boundary<BoundaryType::Periodic, FieldAxes::Z>{};
  }


  template<>
  Boundary<BoundaryType::PML, FieldAxes::X>
  inline createBoundary<BoundaryType::PML, FieldAxes::X, FieldSide::LO>(const Mesh2D& mesh,
                                                                        [[maybe_unused]] size_t nx,
                                                                        [[maybe_unused]] size_t nz,
                                                                        fptype dt,
                                                                        size_t d_pml)
  {
    auto d_sigma = math::linspace(1.0, 0.0, d_pml, false);
    auto d_alpha = math::linspace(0.0, 1.0, d_pml, false);

    auto d_sigma_half = math::meanspace(d_sigma);
    auto d_alpha_half = math::meanspace(d_alpha);

    auto ds_end = d_sigma_half.end() - 1;
    auto da_end = d_alpha_half.end() - 1;
    d_sigma_half.push_back(*(ds_end) - (*(ds_end - 1) - *(ds_end)));
    d_alpha_half.push_back(*(da_end) + (*(da_end) - *(da_end - 1)));

    Boundary<BoundaryType::PML, FieldAxes::X> bc(d_pml, nz);

    auto sigma_bc = calculate_sigma(d_sigma, mesh.dx.front());
    auto kappa_bc = calculate_kappa(d_sigma, kPMLKappaMax[0]);
    auto alpha_bc = calculate_alpha(d_alpha, kPMLAlphaMax[0]);
    calculate_pml_coefs(bc.bE, bc.cE, sigma_bc, kappa_bc, alpha_bc, dt, mesh.dx.front(), d_pml);

    auto sigma_bc_half = calculate_sigma(d_sigma_half, mesh.dx.front());
    auto kappa_bc_half = calculate_kappa(d_sigma_half, kPMLKappaMax[0]);
    auto alpha_bc_half = calculate_alpha(d_alpha_half, kPMLAlphaMax[0]);
    calculate_pml_coefs(bc.bB, bc.cB, sigma_bc_half, kappa_bc_half, alpha_bc_half, dt, mesh.dx.front(), d_pml);

    return bc;
  }

  template<>
  Boundary<BoundaryType::PML, FieldAxes::X>
  inline createBoundary<BoundaryType::PML, FieldAxes::X, FieldSide::HI>(const Mesh2D& mesh,
                                                                        [[maybe_unused]] size_t nx,
                                                                        [[maybe_unused]] size_t nz,
                                                                        fptype dt,
                                                                        size_t d_pml)
  {
    auto d_sigma = math::linspace(1.0, 0.0, d_pml, false);
    auto d_alpha = math::linspace(0.0, 1.0, d_pml, false);
    std::ranges::reverse(d_sigma);
    std::ranges::reverse(d_alpha.begin(), d_alpha.end());

    auto d_sigma_half = math::meanspace(d_sigma);
    auto d_alpha_half = math::meanspace(d_alpha);

    const auto ds_end = d_sigma_half.begin();
    const auto da_end = d_alpha_half.begin();
    d_sigma_half.insert(ds_end, *(ds_end) - (*(ds_end + 1) - *(ds_end)));
    d_alpha_half.insert(da_end, *(da_end) + (*(da_end) - *(da_end + 1)));

    Boundary<BoundaryType::PML, FieldAxes::X> bc(d_pml, nz);

    const auto sigma_bc = calculate_sigma(d_sigma, mesh.dx.back());
    const auto kappa_bc = calculate_kappa(d_sigma, kPMLKappaMax[0]);
    const auto alpha_bc = calculate_alpha(d_alpha, kPMLAlphaMax[0]);
    calculate_pml_coefs(bc.bE, bc.cE, sigma_bc, kappa_bc, alpha_bc, dt, mesh.dx.back(), d_pml);

    const auto sigma_bc_half = calculate_sigma(d_sigma_half, mesh.dx.back());
    const auto kappa_bc_half = calculate_kappa(d_sigma_half, kPMLKappaMax[0]);
    const auto alpha_bc_half = calculate_alpha(d_alpha_half, kPMLAlphaMax[0]);
    calculate_pml_coefs(bc.bB, bc.cB, sigma_bc_half, kappa_bc_half, alpha_bc_half, dt, mesh.dx.back(), d_pml);

    return bc;
  }

  template<>
  Boundary<BoundaryType::PML, FieldAxes::Z>
  inline createBoundary<BoundaryType::PML, FieldAxes::Z, FieldSide::LO>(const Mesh2D& mesh,
                                                                        [[maybe_unused]] size_t nx,
                                                                        [[maybe_unused]] size_t nz,
                                                                        fptype dt,
                                                                        size_t d_pml)
  {
    const auto d_sigma = math::linspace(1.0, 0.0, d_pml, false);
    const auto d_alpha = math::linspace(0.0, 1.0, d_pml, false);

    auto d_sigma_half = math::meanspace(d_sigma);
    auto d_alpha_half = math::meanspace(d_alpha);

    const auto ds_end = d_sigma_half.end() - 1;
    const auto da_end = d_alpha_half.end() - 1;
    d_sigma_half.push_back(*(ds_end) - (*(ds_end - 1) - *(ds_end)));
    d_alpha_half.push_back(*(da_end) + (*(da_end) - *(da_end - 1)));

    Boundary<BoundaryType::PML, FieldAxes::Z> bc(nx, d_pml);

    const auto sigma_bc = calculate_sigma(d_sigma, mesh.dz.front());
    const auto kappa_bc = calculate_kappa(d_sigma, kPMLKappaMax[1]);
    const auto alpha_bc = calculate_alpha(d_alpha, kPMLAlphaMax[1]);
    calculate_pml_coefs(bc.bE, bc.cE, sigma_bc, kappa_bc, alpha_bc, dt, mesh.dz.front(), d_pml);

    const auto sigma_bc_half = calculate_sigma(d_sigma_half, mesh.dz.front());
    const auto kappa_bc_half = calculate_kappa(d_sigma_half, kPMLKappaMax[1]);
    const auto alpha_bc_half = calculate_alpha(d_alpha_half, kPMLAlphaMax[1]);
    calculate_pml_coefs(bc.bB, bc.cB, sigma_bc_half, kappa_bc_half, alpha_bc_half, dt, mesh.dz.front(), d_pml);

    return bc;
  }

  template<>
  Boundary<BoundaryType::PML, FieldAxes::Z>
  inline createBoundary<BoundaryType::PML, FieldAxes::Z, FieldSide::HI>(const Mesh2D& mesh,
                                                                        [[maybe_unused]] size_t nx,
                                                                        [[maybe_unused]] size_t nz,
                                                                        fptype dt,
                                                                        size_t d_pml)
  {
    auto d_sigma = math::linspace(1.0, 0.0, d_pml, false);
    auto d_alpha = math::linspace(0.0, 1.0, d_pml, false);
    std::ranges::reverse(d_sigma.begin(), d_sigma.end());
    std::ranges::reverse(d_alpha.begin(), d_alpha.end());

    auto d_sigma_half = math::meanspace(d_sigma);
    auto d_alpha_half = math::meanspace(d_alpha);

    const auto ds_end = d_sigma_half.begin();
    const auto da_end = d_alpha_half.begin();
    d_sigma_half.insert(ds_end, *(ds_end) - (*(ds_end + 1) - *(ds_end)));
    d_alpha_half.insert(da_end, *(da_end) + (*(da_end) - *(da_end + 1)));

    Boundary<BoundaryType::PML, FieldAxes::Z> bc(nx, d_pml);

    const auto sigma_bc = calculate_sigma(d_sigma, mesh.dz.back());
    const auto kappa_bc = calculate_kappa(d_sigma, kPMLKappaMax[1]);
    const auto alpha_bc = calculate_alpha(d_alpha, kPMLAlphaMax[1]);
    calculate_pml_coefs(bc.bE, bc.cE, sigma_bc, kappa_bc, alpha_bc, dt, mesh.dz.back(), d_pml);

    const auto sigma_bc_half = calculate_sigma(d_sigma_half, mesh.dz.back());
    const auto kappa_bc_half = calculate_kappa(d_sigma_half, kPMLKappaMax[1]);
    const auto alpha_bc_half = calculate_alpha(d_alpha_half, kPMLAlphaMax[1]);
    calculate_pml_coefs(bc.bB, bc.cB, sigma_bc_half, kappa_bc_half, alpha_bc_half, dt, mesh.dz.back(), d_pml);

    return bc;
  }


  /***************************************************/
  /********** Boundary Condition Superclass **********/
  template<BoundaryType X0, BoundaryType X1, BoundaryType Z0, BoundaryType Z1>
  class BoundaryConditions {
  public:
    BoundaryConditions(const Mesh2D& mesh, size_t nx, size_t nz, fptype dt, size_t d_pml = kPMLDepth + params::kNHaloCells)
    : x0(createBoundary<X0, FieldAxes::X, FieldSide::LO>(mesh, nx, nz, dt, d_pml)),
      x1(createBoundary<X1, FieldAxes::X, FieldSide::HI>(mesh, nx, nz, dt, d_pml)),
      z0(createBoundary<Z0, FieldAxes::Z, FieldSide::LO>(mesh, nx, nz, dt, d_pml)),
      z1(createBoundary<Z1, FieldAxes::Z, FieldSide::HI>(mesh, nx, nz, dt, d_pml))
    {}

    Boundary<X0, FieldAxes::X> x0;
    Boundary<X1, FieldAxes::X> x1;
    Boundary<Z0, FieldAxes::Z> z0;
    Boundary<Z1, FieldAxes::Z> z1;
  };


  /*************************************************/
  /********** Boundary Condition Functors **********/
  template<BoundaryType B, FieldAxes F, FieldSide S>
  struct BoundaryFunctor;

  template<FieldSide S>
  struct BoundaryFunctor<BoundaryType::None, FieldAxes::X, S> {
    static void updateE([[maybe_unused]] electromagnetics::emdata::EMData2D<fptype>&, [[maybe_unused]] Boundary<BoundaryType::None, FieldAxes::X>&) {}
    static void updateB([[maybe_unused]] electromagnetics::emdata::EMData2D<fptype>&, Boundary<BoundaryType::None, FieldAxes::X>&, [[maybe_unused]] fptype dt) {}
  };

  template<FieldSide S>
  struct BoundaryFunctor<BoundaryType::None, FieldAxes::Z, S> {
    static void updateE([[maybe_unused]] electromagnetics::emdata::EMData2D<fptype>&, [[maybe_unused]] Boundary<BoundaryType::None, FieldAxes::Z>&) {}
    static void updateB([[maybe_unused]] electromagnetics::emdata::EMData2D<fptype>&, Boundary<BoundaryType::None, FieldAxes::Z>&, [[maybe_unused]] fptype dt) {}
  };

  template<FieldSide S>
  struct BoundaryFunctor<BoundaryType::Periodic, FieldAxes::Z, S> {
    static void updateE(electromagnetics::emdata::EMData2D<fptype>& emdata, [[maybe_unused]] Boundary<BoundaryType::Periodic, FieldAxes::Z>&, [[maybe_unused]] electromagnetics::emdata::Coefficients<fptype>& coefficients)
    {
      auto& [Ex, Ez, By] = emdata.fields;

      const auto& boundary_depth = emdata.boundary_depth;
      const auto num_real_cells = Ex.dims[1] - 2 * params::kNHaloCells;
      const auto lo_bound_idx = params::kNHaloCells;
      const auto hi_bound_idx = (Ex.dims[1] - params::kNHaloCells - 1);

      for (size_t i = 0; i < Ex.dims[0];  i++) {
        Ex(i, params::kNHaloCells) = Ex(i, Ex.dims[1] - params::kNHaloCells - 1);

        for (size_t p = 1; p < params::kNHaloCells - 1; p++) {
          auto pm = p % (num_real_cells - 1);

          auto lo_idx = lo_bound_idx - p;
          Ex(i, lo_idx) = Ex(i, hi_bound_idx - pm);

          auto hi_idx = hi_bound_idx + p;
          Ex(i, hi_idx) = Ex(i, lo_bound_idx + pm);
        }
      }

    }

    static void updateB(electromagnetics::emdata::EMData2D<fptype>& emdata, [[maybe_unused]] Boundary<BoundaryType::Periodic, FieldAxes::Z>&, [[maybe_unused]] fptype dt) {
      auto& [Ex, Ez, By] = emdata.fields;
      for (size_t i = 0; i < By.dims[0]; i++) {
        const auto by_row = By(i, params::kNHaloCells);
        for (size_t k = 0; k < By.dims[1]; k++) {
          By(i, k) = by_row;
        }
      }
    }
  };

  // X0 PML Functor
  template<>
  struct BoundaryFunctor<BoundaryType::PML, FieldAxes::X, FieldSide::LO> {
    static void updateE(electromagnetics::emdata::EMData2D<fptype>& emdata, Boundary<BoundaryType::PML, FieldAxes::X>& x0, electromagnetics::emdata::Coefficients<fptype> coefficients)
    {
      auto& [Ex, Ez, By] = emdata.fields;
      auto& [c_ex_e, c_ex_b, c_ez_e, c_ez_b, c_jx, c_jz, k_x, k_z, k_x_half, k_z_half] = coefficients;
      const auto d_pml = x0.psiE.dims[0];

      // X0
      for(size_t i = 1; i < x0.psiE.dims[0]; ++i) {
        for(size_t k = 0; k < x0.psiE.dims[1]; ++k) {
          auto prev = x0.bE[i] * x0.psiE(i, k);
          auto curl = x0.cE[i] * (By(i, k) - By(i - 1, k));
          x0.psiE(i, k) = prev + curl;
        }
      }
      for(size_t i = 0; i < d_pml; ++i) {
        for(size_t k = 0; k < Ez.dims[1]; ++k) {
          Ez(i, k) += c_ez_b(i, k) * x0.psiE(i, k);
        }
      }
    }

    static void updateB(electromagnetics::emdata::EMData2D<fptype>& emdata, Boundary<BoundaryType::PML, FieldAxes::X>& x0, fptype dt)
    {
      auto& [Ex, Ez, By] = emdata.fields;
      const auto d_pml = x0.psiB.dims[0];

      for (size_t i = 0; i < x0.psiB.dims[0]; ++i) {
        for (size_t k = 0; k < x0.psiB.dims[1]; ++k) {
          const auto prev = x0.bB[i] * x0.psiB(i, k);
          const auto curl = x0.cB[i] * (Ez(i + 1, k) - Ez(i, k));
          x0.psiB(i, k) = prev + curl;
        }
      }
      for (size_t i = 0; i < d_pml; ++i) {
        for (size_t k = 0; k < By.dims[1]; ++k) {
          By(i, k) += dt * x0.psiB(i, k);
        }
      }
    }
  };

  template<>
  struct BoundaryFunctor<BoundaryType::PML, FieldAxes::X, FieldSide::HI> {
    static void updateE(electromagnetics::emdata::EMData2D<fptype>& emdata, Boundary<BoundaryType::PML, FieldAxes::X>& x1, const electromagnetics::emdata::Coefficients<fptype>& coeffs)
    {
      auto& [Ex, Ez, By] = emdata.fields;
      const auto d_pml = x1.psiE.dims[0];

      for(size_t i = 0; i < x1.psiE.dims[0] - 1; ++i) {
        for(size_t k = 0; k < x1.psiE.dims[1]; ++k) {
          size_t offset = i + (By.dims[0] - d_pml);
          auto prev = x1.bE[i] * x1.psiE(i, k);
          auto curl = x1.cE[i] * (By(offset + 1, k) - By(offset, k));
          x1.psiE(i, k) = prev + curl;
        }
      }

      auto start = Ez.dims[0] - d_pml;
      for(size_t i = start; i < Ez.dims[0]; ++i) {
        for(size_t k = 0; k < Ez.dims[1]; ++k) {
          size_t offset = i - (Ez.dims[0] - d_pml);
          Ez(i, k) += coeffs.c_ez_b(i, k) * x1.psiE(offset, k);
        }
      }
    }

    static void updateB(electromagnetics::emdata::EMData2D<fptype>& emdata, Boundary<BoundaryType::PML, FieldAxes::X>& x1, fptype dt)
    {
      auto& [Ex, Ez, By] = emdata.fields;
      const auto d_pml = x1.psiB.dims[0];

      for (size_t i = 0; i < x1.psiB.dims[0]; ++i) {
        for (size_t k = 0; k < x1.psiB.dims[1]; ++k) {
          size_t offset = i + (Ez.dims[0] - d_pml);
          auto prev = x1.bB[i] * x1.psiB(i, k);
          auto curl = x1.cB[i] * (Ez(offset, k) - Ez(offset - 1, k));
          x1.psiB(i, k) = prev + curl;
        }
      }
      auto start = By.dims[0] - d_pml;
      for(size_t i = start; i < By.dims[0]; ++i) {
        for(size_t k = 0; k < By.dims[1]; ++k) {
          size_t offset = i - (By.dims[0] - d_pml);
          By(i, k) += dt * x1.psiB(offset, k);
        }
      }
    }
  };

  template<>
  struct BoundaryFunctor<BoundaryType::PML, FieldAxes::Z, FieldSide::LO> {
    static void updateE(electromagnetics::emdata::EMData2D<fptype>& emdata, Boundary<BoundaryType::PML, FieldAxes::Z>& z0, const electromagnetics::emdata::Coefficients<fptype>& coeffs)
    {
      auto& [Ex, Ez, By] = emdata.fields;
      const auto d_pml = z0.psiE.dims[1];

      for(size_t i = 0; i < z0.psiE.dims[0]; ++i) {
        for(size_t k = 1; k < z0.psiE.dims[1]; ++k) {
          auto prev = z0.bE[k] * z0.psiE(i, k);
          auto curl = z0.cE[k] * (By(i, k) - By(i, k - 1));
          z0.psiE(i, k) = prev + curl;
        }
      }
      for(size_t i = 0; i < Ex.dims[0]; ++i) {
        for(size_t k = 0; k < d_pml; ++k) {
          Ex(i, k) -= coeffs.c_ex_b(i, k) * z0.psiE(i, k);
        }
      }
    }

    static void updateB(electromagnetics::emdata::EMData2D<fptype>& emdata, Boundary<BoundaryType::PML, FieldAxes::Z>& z0, fptype dt)
    {
      auto& [Ex, Ez, By] = emdata.fields;
      const auto d_pml = z0.psiB.dims[1];

      for (size_t i = 0; i < z0.psiB.dims[0]; ++i) {
        for (size_t k = 0; k < z0.psiB.dims[1]; ++k) {
          auto prev = z0.bB[k] * z0.psiB(i, k);
          auto curl = z0.cB[k] * (Ex(i, k + 1) - Ex(i, k));
          z0.psiB(i, k) = prev + curl;
        }
      }
      for(size_t i = 0; i < By.dims[0]; ++i) {
        for(size_t k = 0; k < d_pml; ++k) {
          By(i, k) -= dt * z0.psiB(i, k);
        }
      }
    }
  };

  template<>
  struct BoundaryFunctor<BoundaryType::PML, FieldAxes::Z, FieldSide::HI> {
    static void updateE(electromagnetics::emdata::EMData2D<fptype>& emdata, Boundary<BoundaryType::PML, FieldAxes::Z>& z1, const electromagnetics::emdata::Coefficients<fptype>& coeffs)
    {
      auto& [Ex, Ez, By] = emdata.fields;
      const auto d_pml = z1.psiE.dims[1];

      for(size_t i = 0; i < z1.psiE.dims[0]; ++i) {
        for(size_t k = 0; k < z1.psiE.dims[1] - 1; ++k) {
          const size_t offset = k + (By.dims[1] - d_pml);
          const auto prev = z1.bE[k] * z1.psiE(i, k);
          const auto curl = z1.cE[k] * (By(i, offset + 1) - By(i, offset));
          z1.psiE(i, k) = prev + curl;
        }
      }
      const auto start = Ex.dims[1] - d_pml;
      for(size_t i = 0; i < Ex.dims[0]; ++i) {
        for(size_t k = start; k < Ex.dims[1]; ++k) {
          const size_t offset = k - (Ex.dims[1] - d_pml);
          Ex(i, k) -= coeffs.c_ex_b(i, k) * z1.psiE(i, offset);
        }
      }
    }

    static void updateB(electromagnetics::emdata::EMData2D<fptype>& emdata, Boundary<BoundaryType::PML, FieldAxes::Z>& z1, fptype dt)
    {
      auto& [Ex, Ez, By] = emdata.fields;
      const auto d_pml = z1.psiB.dims[1];

      for (size_t i = 0; i < z1.psiB.dims[0]; ++i) {
        for (size_t k = 0; k < z1.psiB.dims[1]; ++k) {
          const size_t offset = k + (Ex.dims[1] - d_pml);
          const auto prev = z1.bB[k] * z1.psiB(i, k);
          const auto curl = z1.cB[k] * (Ex(i, offset) - Ex(i, offset - 1));
          z1.psiB(i, k) = prev + curl;
        }
      }
      auto start = By.dims[1] - d_pml;
      for(size_t i = 0; i < By.dims[0]; ++i) {
        for(size_t k = start; k < By.dims[1]; ++k) {
          const size_t offset = k - (By.dims[1] - d_pml);
          By(i, k) -= dt * z1.psiB(i, offset);
        }
      }
    }
  };
  //
} // end electromagnetics::boundary_conditions

#endif //TRIFORCE_EM_BOUNDARY_CONDITIONS_H
