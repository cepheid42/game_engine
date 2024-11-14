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

template<typename Derived>
struct Boundary {
  static void updateE(auto&, auto&)      requires is_periodic<Derived> { Derived::updateE(); }
  static void updateH(auto& bc, auto& f) requires is_periodic<Derived> { Derived::updateH(bc, f); }

  static void updateE(auto& bc, auto& f1, auto& f2, const auto& c1)
  requires is_pml<Derived> {
    Derived::updateE(bc, f1, f2, c1);
  }
  static void updateH(auto& bc, auto& f1, auto& f2, const auto& c1)
  requires is_pml<Derived> {
    Derived::updateH(bc, f1, f2, c1);
  }
};

struct Periodic1D : periodic_t<EMFace::X, EMSide::Lo> {
  static void updateE() {}

  static void updateH(auto& bc, auto& f) {
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
struct Periodic2D : periodic_t<F, S> {
  static void updateE() {}

  static void updateH(auto& bc, auto& f) {
    const auto& os = bc.offsets;
    const auto& numInterior = bc.numInterior;
    const auto& hi_idx = bc.hi_idx;
    
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
      }
    }
  }
};

template<EMFace F, EMSide S>
struct Pml2D : pml_t<F, S> {
  static void updateE(auto& bc, auto& f1, const auto& f2, const auto& c1)
  requires (F == EMFace::X)
  {
    const auto y0 = bc.offsets.y0;
    const auto y1 = bc.offsets.y1;

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

      for (size_t j = y0; j < y1; ++j) {
        bc.psi(ipml, j) = bc.b[ipml] * bc.psi(ipml, j) + bc.c[ipml] * (f2(i, j) - f2(i - 1, j));
        f1(i, j) += c1(i, j) * bc.psi(ipml, j);
      }
    }
  }


  static void updateH(auto& bc, auto& f1, const auto& f2, const auto& c1)
  requires (F == EMFace::X)
  {
    const auto y0 = bc.offsets.y0;
    const auto y1 = bc.offsets.y1;

    size_t x0, x1;
    if constexpr (S == EMSide::Lo) {
      x0 = bc.offsets.x0;
      x1 = bc.offsets.x1 - 1;
    } else {
      x0 = bc.offsets.x0 + 1;
      x1 = bc.offsets.x1;
    }

    // DBG(f1.nx(), f1.ny(), x0, x1, y0, y1);

    for (size_t i = x0; i < x1; ++i) {
      size_t ipml;
      if constexpr (S == EMSide::Lo) { ipml = i; }
      else { ipml = i - x0 + 1; }
      // DBG(i, ipml);

      for (size_t j = y0; j < y1; ++j) {
        bc.psi(ipml, j) = bc.b[ipml] * bc.psi(ipml, j) + bc.c[ipml] * (f2(i + 1, j) - f2(i, j));
        f1(i, j) += c1(i, j) * bc.psi(ipml, j);
      }
    }
  }

  static void updateE(auto& bc, auto& f1, const auto& f2, const auto& c1)
  requires (F == EMFace::Y)
  {
    const auto x0 = bc.offsets.x0;
    const auto x1 = bc.offsets.x1;

    size_t y0, y1;
    if constexpr (S == EMSide::Lo) {
      y0 = bc.offsets.y0 + 1;
      y1 = bc.offsets.y1;
    } else {
      y0 = bc.offsets.y0;
      y1 = bc.offsets.y1 - 1;
    }

    for (size_t i = x0; i < x1; ++i) {
      for (size_t j = y0; j < y1; ++j) {
        size_t jpml;
        if constexpr (S == EMSide::Lo) { jpml = j; }
        else { jpml = j - y0; }

        bc.psi(i, jpml) = bc.b[jpml] * bc.psi(i, jpml) + bc.c[jpml] * (f2(i, j) - f2(i, j - 1));
        f1(i, j) -= c1(i, j) * bc.psi(i, jpml);
      }
    }
  }

  static void updateH(auto& bc, auto& f1, const auto& f2, const auto& c1)
  requires (F == EMFace::Y)
  {
    const auto x0 = bc.offsets.x0;
    const auto x1 = bc.offsets.x1;

    size_t y0, y1;
    if constexpr (S == EMSide::Lo) {
      y0 = bc.offsets.y0;
      y1 = bc.offsets.y1 - 1;
    } else {
      y0 = bc.offsets.y0 + 1;
      y1 = bc.offsets.y1;
    }

    // DBG(f1.nx(), f1.ny(), x0, x1, y0, y1);

    for (size_t i = x0; i < x1; ++i) {
      for (size_t j = y0; j < y1; ++j) {
        size_t jpml;
        if constexpr (S == EMSide::Lo) { jpml = j; }
        else { jpml = j - y0 + 1; }

        // DBG(j, jpml);

        bc.psi(i, jpml) = bc.b[jpml] * bc.psi(i, jpml) + bc.c[jpml] * (f2(i, j + 1) - f2(i, j));
        f1(i, j) -= c1(i, j) * bc.psi(i, jpml);
      }
    }
  }
};




#endif //BOUNDARIES_H
