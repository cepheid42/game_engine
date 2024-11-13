//
// Created by cepheid on 10/17/24.
//

#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include "electromagnetics.param"

template<EMFace F>
struct periodic_t {
  static constexpr EMFace face = F;
};

template<EMFace F>
struct pml_t {
  static constexpr EMFace face = F;
};

template<typename T>
concept is_periodic = std::derived_from<T, periodic_t<T::face>>;

template<typename T>
concept is_pml = std::derived_from<T, pml_t<T::face>>;


template<typename Derived>
struct Boundary {
  static void updateE(const auto&, const auto&) requires is_periodic<Derived> { Derived::updateE(); }
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


struct Periodic1D : periodic_t<EMFace::X> {
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

template<EMFace F>
struct Periodic2D : periodic_t<F> {
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
//
// struct Periodic3D : Boundary<Periodic3D> {
//   void updateE() {}
//   void updateH() {}
// };


// struct Pml1D : pml_bc {
//   static void updateE() {}
//
//   static void updateH() {}
// };

template<EMFace F>
struct Pml2D : pml_t<F> {
  static void updateE(auto& bc, auto& f1, const auto& f2, const auto& c1)
  requires (F == EMFace::X)
  {
    const auto& os = bc.offsets;

    for (size_t i = os.x0 + 1; i < os.x1; ++i) {
      const auto ipml = i - (os.x0 + 1);
      for (size_t j = os.y0; j < os.y1; ++j) {
        bc.psi(ipml, j) = bc.b[ipml] * bc.psi(ipml, j) + bc.c[ipml] * (f2(i, j) - f2(i - 1, j));
        f1(i, j) += c1(i, j) * bc.psi(ipml, j);
      }
    }
  }

  static void updateH(auto& bc, auto& f1, const auto& f2, const auto& c1)
  requires (F == EMFace::X)
  {
    const auto& os = bc.offsets;

    for (size_t i = os.x0; i < os.x1; ++i) {
      const auto ipml = i - os.x0;
      for (size_t j = os.y0; j < os.y1; ++j) {
        bc.psi(ipml, j) = bc.b[ipml] * bc.psi(ipml, j) + bc.c[ipml] * (f2(i + 1, j) - f2(i, j));
        f1(i, j) -= c1(i, j) * bc.psi(ipml, j);
      }
    }
  }
};

// struct Pml3D : Boundary<Pml3D> {
//   void updateE() {}
//   void updateH() {}
// };
//
//



#endif //BOUNDARIES_H
