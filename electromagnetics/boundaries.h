//
// Created by cepheid on 10/17/24.
//

#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include "electromagnetics.param"

template<EMFace F, EMSide S>
struct periodic_t {
  static constexpr EMFace face = F;
  static constexpr EMSide side = S;
};

template<EMFace F, EMSide S>
struct pml_t {
  static constexpr EMFace face = F;
  static constexpr EMSide side = S;
};

template<typename T>
concept is_periodic = std::derived_from<T, periodic_t<T::face, T::side>>;

template<typename T>
concept is_pml = std::derived_from<T, pml_t<T::face, T::side>>;

template<typename UpdateFunc>
struct Boundary {
  static void updateE(auto& bc, auto& f1, auto& f2, const auto& c1) {
    if constexpr (is_pml<UpdateFunc>) {
      // DBG("PmlE()");
      UpdateFunc::updateE(bc, f1, f2, c1);
    }
    // Periodic & Reflecting do not update the E-field components
  }

  static void updateH(auto& bc, auto& f1, auto& f2, const auto& c1) {
    if constexpr (is_periodic<UpdateFunc>) {
      // DBG("PeriodicH()");
      UpdateFunc::updateH(bc, f1);
    } else if constexpr (is_pml<UpdateFunc>) {
      // DBG("PmlH()");
      UpdateFunc::updateH(bc, f1, f2, c1);
    }
    // Reflecting does not update H field components
  }
};

struct ReflectingBCUpdate {};

template<EMFace F, EMSide S>
struct Periodic1DUpdate : periodic_t<F, S> {
  // static void updateE() {}
  static void updateH(auto&, auto&) { DBG("Periodic::updateH(hi)"); }
  static void updateH(auto& bc, auto& f) requires (S == EMSide::Lo) {
    DBG("Periodic::updateH(lo)");
    const auto& os = bc.offsets;
    const auto& numInterior = bc.numInterior;
    const auto& hi_idx = bc.hi_idx;

    for (size_t i = os.x0; i < os.x1; ++i) {
      const auto pm = i % numInterior;

      f[bc.depth - 1 - i] = f[hi_idx - pm];
      f[hi_idx + 1 + i] = f[bc.depth + pm];
    }
  }
};

template<EMFace F, EMSide S>
struct Periodic2DUpdate : periodic_t<F, S> {
  // static void updateE() {}
  static void updateH(auto&, auto&) {} // Only lo-side periodic is used
  static void updateH(auto& bc, auto& f) requires (S == EMSide::Lo) {
    const auto& os = bc.offsets;
    const auto& numInterior = bc.numInterior;
    const auto& hi_idx = bc.hi_idx;

    // DBG("Periodic2D::updateH");

    for (size_t i = os.x0; i < os.x1; ++i) {
      for (size_t j = os.y0; j < os.y1; ++j) {
        if constexpr (F == EMFace::X) {
          const auto pm = i % numInterior;

          f(bc.depth - 1 - i, j) = f(hi_idx - pm, j);
          f(hi_idx + 1 + i, j) = f(bc.depth + pm, j);
        } else {
          const auto pm = j % numInterior;

          f(i, bc.depth - 1 - j) = f(i, hi_idx - pm);
          f(i, hi_idx + 1 + j) = f(i, bc.depth + pm);
        }
      } // end j-loop
    } // end i-loop
  } // end updateH
};

template<EMFace F, EMSide S>
struct Periodic3DUpdate : periodic_t<F, S> {
  // static void updateE() { DBG("Periodic3D::updateE()"); }
  static void updateH(auto&, auto&) {} // Only lo-side periodic is used
  static void updateH(auto& bc, auto& f) requires (S == EMSide::Lo) {
    const auto& os = bc.offsets;
    const auto& numInterior = bc.numInterior;
    const auto& hi_idx = bc.hi_idx;

    for (size_t i = os.x0; i < os.x1; ++i) {
      for (size_t j = os.y0; j < os.y1; ++j) {
        for (size_t k = os.z0; k < os.z1; ++k) {
          if constexpr (F == EMFace::X) {
            const auto pm = i % numInterior;

            f(bc.depth - 1 - i, j, k) = f(hi_idx - pm, j, k);
            f(hi_idx + 1 + i, j, k) = f(bc.depth + pm, j, k);
          } else if constexpr (F == EMFace::Y) {
            const auto pm = j % numInterior;

            f(i, bc.depth - 1 - j, k) = f(i, hi_idx - pm, k);
            f(i, hi_idx + 1 + j, k) = f(i, bc.depth + pm, k);
          } else {
            const auto pm = k % numInterior;

            f(i, j, bc.depth - 1 - k) = f(i, j, hi_idx - pm);
            f(i, j, hi_idx + 1 + k) = f(i, j, bc.depth + pm);
          }
        } // end k-loop
      } // end j-loop
    } // end i-loop
  } // end updateH
};

template<EMFace F, EMSide S>
struct Pml1DUpdate : pml_t<F, S> {
  static void updateE(auto& bc, auto& f1, const auto& f2, const auto& c1) {
    size_t x0, x1;
    if constexpr (S == EMSide::Lo) {
      x0 = bc.offsets.x0 + 1;
      x1 = bc.offsets.x1;
    } else {
      x0 = bc.offsets.x0;
      x1 = bc.offsets.x1 - 1;
    }

    for (size_t i = x0; i < x1; ++i) {
      size_t ipml;
      if constexpr (S == EMSide::Lo) { ipml = i; }
      else { ipml = i - x0; }

        bc.psi[ipml] = bc.b[ipml] * bc.psi[ipml] + bc.c[ipml] * (f2[i] - f2[i - 1]);
        f1[i] += c1[i] * bc.psi[ipml];
    } // end i-loop
  } // end updateE

  static void updateH(auto& bc, auto& f1, const auto& f2, const auto& c1) {
    size_t x0, x1;
    if constexpr (S == EMSide::Lo) {
      x0 = bc.offsets.x0;
      x1 = bc.offsets.x1 - 1;
    } else {
      x0 = bc.offsets.x0 + 1;
      x1 = bc.offsets.x1;
    }

    for (size_t i = x0; i < x1; ++i) {
      size_t ipml;
      if constexpr (S == EMSide::Lo) { ipml = i; }
      else { ipml = i - x0 + 1; }

      bc.psi[ipml] = bc.b[ipml] * bc.psi[ipml] + bc.c[ipml] * (f2[i + 1] - f2[i]);
      f1[i] += c1[i] * bc.psi[ipml];
    } // end i-loop
  } // end updateH
};

template<EMFace F, EMSide S, bool Negate>
struct Pml2DUpdate : pml_t<F, S> {
  static void updateE(auto& bc, auto& f1, const auto& f2, const auto& c1)
  requires (F == EMFace::X)
  {
    const auto& y0 = bc.offsets.y0;
    const auto& y1 = bc.offsets.y1;

    size_t x0, x1;
    if constexpr (S == EMSide::Lo) {
      // DBG("Pml2DX::updateE(lo)");
      x0 = bc.offsets.x0 + 1;
      x1 = bc.offsets.x1;
    } else {
      // DBG("Pml2DX::updateE(hi)");
      x0 = bc.offsets.x0;
      x1 = bc.offsets.x1 - 1;
    }

    size_t ipml;
#pragma omp parallel for collapse(2) num_threads(NTHREADS_BC) default(shared) private(ipml)
    for (size_t i = x0; i < x1; ++i) {
      for (size_t j = y0; j < y1; ++j) {
        if constexpr (S == EMSide::Lo) { ipml = i; }
        else { ipml = i - x0; }

        bc.psi(ipml, j) = bc.b[ipml] * bc.psi(ipml, j) + bc.c[ipml] * (f2(i, j) - f2(i - 1, j));
        if constexpr (Negate) {
          f1(i, j) -= c1(i, j) * bc.psi(ipml, j);
        } else {
          f1(i, j) += c1(i, j) * bc.psi(ipml, j);
        }
      } // end j-loop
    } // end i-loop
  } // end updateE (X-Face)

  static void updateH(auto& bc, auto& f1, const auto& f2, const auto& c1)
  requires (F == EMFace::X)
  {
    const auto& y0 = bc.offsets.y0;
    const auto& y1 = bc.offsets.y1;

    size_t x0, x1;
    if constexpr (S == EMSide::Lo) {
      // DBG("Pml2DX::updateH(lo)");
      x0 = bc.offsets.x0;
      x1 = bc.offsets.x1 - 1;
    } else {
      // DBG("Pml2DX::updateH(hi)");
      x0 = bc.offsets.x0 + 1;
      x1 = bc.offsets.x1;
    }
    size_t ipml;
#pragma omp parallel for collapse(2) num_threads(NTHREADS_BC) default(shared) private(ipml)
    for (size_t i = x0; i < x1; ++i) {
      for (size_t j = y0; j < y1; ++j) {
        if constexpr (S == EMSide::Lo) { ipml = i; }
        else { ipml = i - x0 + 1; }

        bc.psi(ipml, j) = bc.b[ipml] * bc.psi(ipml, j) + bc.c[ipml] * (f2(i + 1, j) - f2(i, j));
        if constexpr (Negate) {
          f1(i, j) -= c1(i, j) * bc.psi(ipml, j);
        } else {
          f1(i, j) += c1(i, j) * bc.psi(ipml, j);
        }
      } // end j-loop
    } // end i-loop
  } // end updateH (X-Face)

  static void updateE(auto& bc, auto& f1, const auto& f2, const auto& c1)
  requires (F == EMFace::Y)
  {
    const auto& x0 = bc.offsets.x0;
    const auto& x1 = bc.offsets.x1;

    size_t y0, y1;
    if constexpr (S == EMSide::Lo) {
      // DBG("Pml2DY::updateE(lo)");
      y0 = bc.offsets.y0 + 1;
      y1 = bc.offsets.y1;
    } else {
      // DBG("Pml2DY::updateE(hi)");
      y0 = bc.offsets.y0;
      y1 = bc.offsets.y1 - 1;
    }

    size_t jpml;
#pragma omp parallel for collapse(2) num_threads(NTHREADS_BC) default(shared) private(jpml)
    for (size_t i = x0; i < x1; ++i) {
      for (size_t j = y0; j < y1; ++j) {
        if constexpr (S == EMSide::Lo) { jpml = j; }
        else { jpml = j - y0; }

        bc.psi(i, jpml) = bc.b[jpml] * bc.psi(i, jpml) + bc.c[jpml] * (f2(i, j) - f2(i, j - 1));
        if constexpr (Negate) {
          f1(i, j) -= c1(i, j) * bc.psi(i, jpml);
        } else {
          f1(i, j) += c1(i, j) * bc.psi(i, jpml);
        }
      } // end j-loop
    } // end i-loop
  } // end updateE (Y-Face)

  static void updateH(auto& bc, auto& f1, const auto& f2, const auto& c1)
  requires (F == EMFace::Y)
  {
    const auto& x0 = bc.offsets.x0;
    const auto& x1 = bc.offsets.x1;

    size_t y0, y1;
    if constexpr (S == EMSide::Lo) {
      // DBG("Pml2DY::updateH(lo)");
      y0 = bc.offsets.y0;
      y1 = bc.offsets.y1 - 1;
    } else {
      // DBG("Pml2DY::updateH(hi)");
      y0 = bc.offsets.y0 + 1;
      y1 = bc.offsets.y1;
    }

    size_t jpml;
#pragma omp parallel for collapse(2) num_threads(NTHREADS_BC) default(shared) private(jpml)
    for (size_t i = x0; i < x1; ++i) {
      for (size_t j = y0; j < y1; ++j) {
        if constexpr (S == EMSide::Lo) { jpml = j; }
        else { jpml = j - y0 + 1; }

        bc.psi(i, jpml) = bc.b[jpml] * bc.psi(i, jpml) + bc.c[jpml] * (f2(i, j + 1) - f2(i, j));
        if constexpr (Negate) {
          f1(i, j) -= c1(i, j) * bc.psi(i, jpml);
        } else {
          f1(i, j) += c1(i, j) * bc.psi(i, jpml);
        }
      } // end j-loop
    } // end i-loop
  } // end updateH (Y-Face)
};

template<EMFace F, EMSide S, bool Negate>
struct Pml3DUpdate : pml_t<F, S> {
  static void updateE(auto& bc, auto& f1, const auto& f2, const auto& c1)
  requires (F == EMFace::X)
  {
    const auto& y0 = bc.offsets.y0;
    const auto& y1 = bc.offsets.y1;
    const auto& z0 = bc.offsets.z0;
    const auto& z1 = bc.offsets.z1;

    size_t x0, x1;
    if constexpr (S == EMSide::Lo) {
      x0 = bc.offsets.x0 + 1;
      x1 = bc.offsets.x1;
    } else {
      x0 = bc.offsets.x0;
      x1 = bc.offsets.x1 - 1;
    }

    size_t ipml;
#pragma omp parallel for collapse(3) num_threads(NTHREADS_BC) default(shared) private(ipml)
    for (size_t i = x0; i < x1; ++i) {
      for (size_t j = y0; j < y1; ++j) {
        for (size_t k = z0; k < z1; ++k) {
          if constexpr (S == EMSide::Lo) { ipml = i; }
          else { ipml = i - x0; }

          bc.psi(ipml, j, k) = bc.b[ipml] * bc.psi(ipml, j, k) + bc.c[ipml] * (f2(i, j, k) - f2(i - 1, j, k));
          if constexpr (Negate) {
            f1(i, j, k) -= c1(i, j, k) * bc.psi(ipml, j, k);
          } else {
            f1(i, j, k) += c1(i, j, k) * bc.psi(ipml, j, k);
          }
        } // end k-loop
      } // end j-loop
    } // end i-loop
  } // end updateE (X-Face)

  static void updateH(auto& bc, auto& f1, const auto& f2, const auto& c1)
  requires (F == EMFace::X)
  {
    const auto& y0 = bc.offsets.y0;
    const auto& y1 = bc.offsets.y1;
    const auto& z0 = bc.offsets.z0;
    const auto& z1 = bc.offsets.z1;

    size_t x0, x1;
    if constexpr (S == EMSide::Lo) {
      x0 = bc.offsets.x0;
      x1 = bc.offsets.x1 - 1;
    } else {
      x0 = bc.offsets.x0 + 1;
      x1 = bc.offsets.x1;
    }

    size_t ipml;
#pragma omp parallel for collapse(3) num_threads(NTHREADS_BC) default(shared) private(ipml)
    for (size_t i = x0; i < x1; ++i) {
      for (size_t j = y0; j < y1; ++j) {
        for (size_t k = z0; k < z1; ++k) {
          if constexpr (S == EMSide::Lo) { ipml = i; }
          else { ipml = i - x0 + 1; }

          bc.psi(ipml, j, k) = bc.b[ipml] * bc.psi(ipml, j, k) + bc.c[ipml] * (f2(i + 1, j, k) - f2(i, j, k));
          if constexpr (Negate) {
            f1(i, j, k) -= c1(i, j, k) * bc.psi(ipml, j, k);
          } else {
            f1(i, j, k) += c1(i, j, k) * bc.psi(ipml, j, k);
          }
        } // end k-loop
      } // end j-loop
    } // end i-loop
  } // end updateH (X-Face)

  static void updateE(auto& bc, auto& f1, const auto& f2, const auto& c1)
  requires (F == EMFace::Y)
  {
    const auto& x0 = bc.offsets.x0;
    const auto& x1 = bc.offsets.x1;
    const auto& z0 = bc.offsets.z0;
    const auto& z1 = bc.offsets.z1;

    size_t y0, y1;
    if constexpr (S == EMSide::Lo) {
      y0 = bc.offsets.y0 + 1;
      y1 = bc.offsets.y1;
    } else {
      y0 = bc.offsets.y0;
      y1 = bc.offsets.y1 - 1;
    }

    size_t jpml;
#pragma omp parallel for collapse(3) num_threads(NTHREADS_BC) default(shared) private(jpml)
    for (size_t i = x0; i < x1; ++i) {
      for (size_t j = y0; j < y1; ++j) {
        for (size_t k = z0; k < z1; ++k) {
          if constexpr (S == EMSide::Lo) { jpml = j; }
          else { jpml = j - y0; }

          bc.psi(i, jpml, k) = bc.b[jpml] * bc.psi(i, jpml, k) + bc.c[jpml] * (f2(i, j, k) - f2(i, j - 1, k));
          if constexpr (Negate) {
            f1(i, j, k) -= c1(i, j, k) * bc.psi(i, jpml, k);
          } else {
            f1(i, j, k) += c1(i, j, k) * bc.psi(i, jpml, k);
          }
        } // end k-loop
      } // end j-loop
    } // end i-loop
  } // end updateE (Y-Face)

  static void updateH(auto& bc, auto& f1, const auto& f2, const auto& c1)
  requires (F == EMFace::Y)
  {
    const auto& x0 = bc.offsets.x0;
    const auto& x1 = bc.offsets.x1;
    const auto& z0 = bc.offsets.z0;
    const auto& z1 = bc.offsets.z1;

    size_t y0, y1;
    if constexpr (S == EMSide::Lo) {
      y0 = bc.offsets.y0;
      y1 = bc.offsets.y1 - 1;
    } else {
      y0 = bc.offsets.y0 + 1;
      y1 = bc.offsets.y1;
    }

    size_t jpml;
#pragma omp parallel for collapse(3) num_threads(NTHREADS_BC) default(shared) private(jpml)
    for (size_t i = x0; i < x1; ++i) {
      for (size_t j = y0; j < y1; ++j) {
        for (size_t k = z0; k < z1; ++k) {
          if constexpr (S == EMSide::Lo) { jpml = j; }
          else { jpml = j - y0 + 1; }

          bc.psi(i, jpml, k) = bc.b[jpml] * bc.psi(i, jpml, k) + bc.c[jpml] * (f2(i, j + 1, k) - f2(i, j, k));
          if constexpr (Negate) {
            f1(i, j, k) -= c1(i, j, k) * bc.psi(i, jpml, k);
          } else {
            f1(i, j, k) += c1(i, j, k) * bc.psi(i, jpml, k);
          }
        } // end k-loop
      } // end j-loop
    } // end i-loop
  } // end updateH (Y-Face)

  static void updateE(auto& bc, auto& f1, const auto& f2, const auto& c1)
  requires (F == EMFace::Z)
  {
    const auto& x0 = bc.offsets.x0;
    const auto& x1 = bc.offsets.x1;
    const auto& y0 = bc.offsets.y0;
    const auto& y1 = bc.offsets.y1;

    size_t z0, z1;
    if constexpr (S == EMSide::Lo) {
      z0 = bc.offsets.z0 + 1;
      z1 = bc.offsets.z1;
    } else {
      z0 = bc.offsets.z0;
      z1 = bc.offsets.z1 - 1;
    }

    size_t kpml;
#pragma omp parallel for collapse(3) num_threads(16) default(shared) private(kpml)
    for (size_t i = x0; i < x1; ++i) {
      for (size_t j = y0; j < y1; ++j) {
        for (size_t k = z0; k < z1; ++k) {
          if constexpr (S == EMSide::Lo) { kpml = k; }
          else { kpml = k - z0; }

          bc.psi(i, j, kpml) = bc.b[kpml] * bc.psi(i, j, kpml) + bc.c[kpml] * (f2(i, j, k) - f2(i, j, k - 1));
          if constexpr (Negate) {
            f1(i, j, k) -= c1(i, j, k) * bc.psi(i, j, kpml);
          } else {
            f1(i, j, k) += c1(i, j, k) * bc.psi(i, j, kpml);
          }
        } // end k-loop
      } // end j-loop
    } // end i-loop
  } // end updateE (Z-Face)
  
  static void updateH(auto& bc, auto& f1, const auto& f2, const auto& c1)
  requires (F == EMFace::Z)
  {
    const auto& x0 = bc.offsets.x0;
    const auto& x1 = bc.offsets.x1;
    const auto& y0 = bc.offsets.y0;
    const auto& y1 = bc.offsets.y1;

    size_t z0, z1;
    if constexpr (S == EMSide::Lo) {
      z0 = bc.offsets.z0;
      z1 = bc.offsets.z1 - 1;
    } else {
      z0 = bc.offsets.z0 + 1;
      z1 = bc.offsets.z1;
    }

    size_t kpml;
#pragma omp parallel for collapse(3) num_threads(16) default(shared) private(kpml)
    for (size_t i = x0; i < x1; ++i) {
      for (size_t j = y0; j < y1; ++j) {
        for (size_t k = z0; k < z1; ++k) {
          if constexpr (S == EMSide::Lo) { kpml = k; }
          else { kpml = k - z0 + 1; }

          bc.psi(i, j, kpml) = bc.b[kpml] * bc.psi(i, j, kpml) + bc.c[kpml] * (f2(i, j, k + 1) - f2(i, j, k));
          if constexpr (Negate) {
            f1(i, j, k) -= c1(i, j, k) * bc.psi(i, j, kpml);
          } else {
            f1(i, j, k) += c1(i, j, k) * bc.psi(i, j, kpml);
          }
        } // end k-loop
      } // end j-loop
    } // end i-loop
  } // end updateH (Z-Face)
};

#endif //BOUNDARIES_H
