//
// Created by cepheid on 9/9/24.
//

#ifndef OLD_FUNCTORS_H
#define OLD_FUNCTORS_H

#include <cstddef>

//using fptype = double;
//
//template<size_t... Na>
//struct Array {
//  fptype data[(Na * ...)]{};
//  size_t shape[sizeof...(Na)] = {Na...};
//
//  [[nodiscard]] auto get_index(size_t i, size_t j) const { return j + (shape[1] * i); }
//  [[nodiscard]] auto get_index(size_t i, size_t j, size_t k) const { return k + (shape[1] * j) + (shape[2] * shape[1] * i); }
//
//  // 1D indexing
//  auto& operator()(const size_t i) { return data[i]; }
//  const auto& operator()(const size_t i) const { return data[i]; }
//
//  // 2D indexing
//  auto& operator()(const size_t i, const size_t j) { assert(i < shape[0] and j < shape[1]); return data[get_index(i, j)]; }
//  const auto& operator()(const size_t i, const size_t j) const { assert(i < shape[0] and j < shape[1]); return data[get_index(i, j)]; }
//
//  // 3D indexing
//  auto& operator()(const size_t i, const size_t j, const size_t k) { return data[get_index(i, j, k)]; }
//  const auto& operator()(const size_t i, const size_t j, const size_t k) const { return data[get_index(i, j, k)]; }
//
//  // [[nodiscard]] auto diff_x(size_t i) const { return data[i + 1] - data[i]; }
//  [[nodiscard]] auto diff_x(size_t i, size_t j) const { return data[get_index(i, j)] - data[get_index(i - 1, j)]; }
//  // [[nodiscard]] auto diff_x(size_t i, size_t j, size_t k) const { return data[get_index(i + 1, j, k)] - data[get_index(i, j, k)]; }
//
//  // [[nodiscard]] auto diff_y(size_t i) const { return data[i + 1] - data[i]; }
//  [[nodiscard]] auto diff_y(size_t i, size_t j) const { return data[get_index(i, j + 1)] - data[get_index(i, j)]; }
//  // [[nodiscard]] auto diff_y(size_t i, size_t j, size_t k) const { return data[get_index(i, j + 1, k)] - data[get_index(i, j, k)]; }
//
//  // [[nodiscard]] auto diff_z(size_t i) const { return data[i + 1] - data[i]; }
//  // [[nodiscard]] auto diff_z(size_t i, size_t j) const { return data[get_index(i, j + 1)] - data[get_index(i, j)]; }
//  // [[nodiscard]] auto diff_z(size_t i, size_t j, size_t k) const { return data[get_index(i, j, k + 1)] - data[get_index(i, j, k)]; }
//};
//
enum class Derivative{ DX, DY, DZ, NoOp };

struct EMPTYARRAY {
  constexpr auto operator()(size_t...) const { return 0.0; }
};

template<size_t DIM>
struct IntegratorOffsets {
  std::array<size_t, 2 * DIM> offsets;
};

// inline constexpr IntegratorOffsets<1> Ex1D_Offsets{1u, 1u};
// inline constexpr IntegratorOffsets<1> Ey1D_Offsets = Ex1D_Offsets;
// inline constexpr IntegratorOffsets<1> Ez1D_Offsets = Ex1D_Offsets;

// inline constexpr IntegratorOffsets<2> Ex2D_Offsets{1u, 1u, 1u, 1u};
// inline constexpr IntegratorOffsets<2> Ey2D_Offsets{0u, 0u, };
// inline constexpr IntegratorOffsets<2> Ez2D_Offsets{0u, 0u};

//====== Curl Operators =======
//=============================
template<Derivative D, typename... IDXS>
struct curl {
  static constexpr auto apply(const auto&, IDXS...) { return 0.0; }
};

template<typename... IDXS>
struct curl<Derivative::DX, IDXS...> {
  static auto apply(const auto& f, IDXS... idxs) {
    return f.diff_x(idxs...);
  }
};

template<typename... IDXS>
struct curl<Derivative::DY, IDXS...> {
  static auto apply(const auto& f, IDXS... idxs) { return f.diff_y(idxs...); }
};

template<typename... IDXS>
struct curl<Derivative::DZ, IDXS...> {
  static auto apply(const auto& f, IDXS... idxs) {
    return f.diff_z(idxs...);
  }
};

//=================== Field Functors ========================
//===========================================================
template<typename TL, Derivative ACURL, Derivative BCURL, typename... IDXS>
struct UpdateFunctor {
  using curl1 = curl<ACURL, IDXS...>;
  using curl2 = curl<BCURL, IDXS...>;
  using F = TypeListAt<0, TL>;
  using D1 = TypeListAt<1, TL>;
  using D2 = TypeListAt<2, TL>;
  using S = TypeListAt<3, TL>;
  using C1 = TypeListAt<4, TL>;
  using C2 = TypeListAt<5, TL>;
  using C3 = TypeListAt<6, TL>;

  static auto apply(F& f, const D1& d1, const D2& d2, const S& src, const C1& c_f, const C2& c_d, const C3& c_src, IDXS... idxs) {
    const auto self = c_f(idxs...) * f(idxs...);
    const auto diff1 = curl1::apply(d1, idxs...);
    const auto diff2 = curl2::apply(d2, idxs...);
    const auto diff = c_d(idxs...) * (diff1 - diff2);
    const auto current = c_src(idxs...) * src(idxs...);
    DBG(self, diff1, diff2, current);
    return self + diff - current;
  }
};



template<typename TL, Derivative ACURL, Derivative BCURL>
struct FieldIntegrator1D {
  using F = TypeListAt<0, TL>;
  using D1 = TypeListAt<1, TL>;
  using D2 = TypeListAt<2, TL>;
  using J1 = TypeListAt<3, TL>;
  using C1 = TypeListAt<4, TL>;
  using C2 = TypeListAt<5, TL>;
  using C3 = TypeListAt<6, TL>;
  using C4 = TypeListAt<7, TL>;
  using update_func = UpdateFunctor<TL, ACURL, BCURL, std::size_t>;
  static auto apply(F& f, const D1& d1, const D2& d2, const J1& j1, const C1& c1, const C2& c2, const C3& c3, const C4& c4) {
    for (size_t i = 0; i < f.shape[0]; i++) {
      update_func::apply(f, d1, d2, j1, c1, c2, c3, c4, i);
    }
  }
};

template<typename TL, Derivative ACURL, Derivative BCURL>
struct FieldIntegrator2D {
  using F = TypeListAt<0, TL>;
  using D1 = TypeListAt<1, TL>;
  using D2 = TypeListAt<2, TL>;
  using S = TypeListAt<3, TL>;
  using C1 = TypeListAt<4, TL>;
  using C2 = TypeListAt<5, TL>;
  using C3 = TypeListAt<6, TL>;
  using update_func = UpdateFunctor<TL, ACURL, BCURL, std::size_t, std::size_t>;

  static auto apply(F& f, const D1& d1, const D2& d2, const S& src, const C1& c_f, const C2& c_d, const C3& c_src, const IntegratorOffsets<2>& o) {
    const auto [x0, x1, y0, y1] = o.offsets;

    for (size_t i = x0; i < f.shape[0] - x1; i++) {
      for (size_t j = y0; j < f.shape[1] - y1; j++) {
        DBG(i, j, f.shape[0], f.shape[1]);
        update_func::apply(f, d1, d2, src, c_f, c_d, c_src, i, j);
      }
    }
  }
};

template<typename TL, Derivative ACURL, Derivative BCURL>
struct FieldIntegrator3D {
  using F = TypeListAt<0, TL>;
  using D1 = TypeListAt<1, TL>;
  using D2 = TypeListAt<2, TL>;
  using J1 = TypeListAt<3, TL>;
  using C1 = TypeListAt<4, TL>;
  using C2 = TypeListAt<5, TL>;
  using C3 = TypeListAt<6, TL>;
  using C4 = TypeListAt<7, TL>;
  using update_func = UpdateFunctor<TL, ACURL, BCURL, std::size_t, std::size_t, std::size_t>;
  static auto apply(F& f, const D1& d1, const D2& d2, const J1& j1, const C1& c1, const C2& c2, const C3& c3, const C4& c4) {
    for (size_t i = 0; i < f.shape[0]; i++) {
      for (size_t j = 0; j < f.shape[1]; j++) {
        for (size_t k = 0; k < f.shape[2]; k++) {
          update_func::apply(f, d1, d2, j1, c1, c2, c3, c4, i, j, k);
        }
      }
    }
  }
};

struct FieldIntegratorNull { static auto apply(...) {} };


template<typename EIX, typename EIY, typename EIZ, typename BIX, typename BIY, typename BIZ>
struct Electromagnetics {
  static void updateE(auto& fields, const auto& coeffs) {
    constexpr IntegratorOffsets<2> zero_offsets = {0u, 0u, 0u, 0u};
    constexpr IntegratorOffsets<2> one_offsets = {1u, 1u, 1u, 1u};

    // EIX::apply(fields.Ex, fields.Bz, fields.By, fields.Jx, coeffs.cexe, coeffs.cexh, coeffs.cjx, zero_offsets);
    // EIY::apply(fields.Ey, fields.Bx, fields.Bz, fields.Jy, coeffs.ceye, coeffs.ceyh, coeffs.cjy, zero_offsets);
    EIZ::apply(fields.Ez, fields.By, fields.Bx, fields.Jz, coeffs.ceze, coeffs.cezh, coeffs.cjz, one_offsets);
  }

  static void updateB(auto& fields, const auto& coeffs) {
    constexpr IntegratorOffsets<2> zero_offsets = {0u, 0u, 0u, 0u};
    BIX::apply(fields.Bx, fields.Ey, fields.Ez, EMPTYARRAY{}, coeffs.chxh, coeffs.chxe, EMPTYARRAY{}, zero_offsets);
    BIY::apply(fields.By, fields.Ez, fields.Ex, EMPTYARRAY{}, coeffs.chyh, coeffs.chye, EMPTYARRAY{}, zero_offsets);
    BIZ::apply(fields.Bz, fields.Ex, fields.Ey, EMPTYARRAY{}, coeffs.chzh, coeffs.chze, EMPTYARRAY{}, zero_offsets);
  }
};


#endif //OLD_FUNCTORS_H
