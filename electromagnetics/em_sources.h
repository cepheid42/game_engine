//
// Created by cepheid on 12/13/24.
//

#ifndef EM_SOURCES_H
#define EM_SOURCES_H

#include <cmath>

#include "electromagnetics.param"

namespace tf::electromagnetics::sources
{
  template<typename T>
  struct TemporalSource {
    virtual ~TemporalSource() = default;
    [[nodiscard]] virtual T eval(T) const = 0;
  };

  template<typename T>
  struct GaussianSource final : TemporalSource<T> {
    [[nodiscard]] T eval(const T t) const override {
      constexpr auto tol = 1e-12;
      const auto val = std::exp(-1.0 * std::pow((t - delay) / width, power));
      return val > tol ? val : 0.0;
    }

    T width;
    T power;
    // T cutoff;
    T delay;
  };

  template<typename T>
  struct RickerSource final : TemporalSource<T> {
    [[nodiscard]] T eval(const T q) const override {
      constexpr auto Np = 20.0;
      constexpr auto Md = 2.0;

      const auto alpha = (M_PI * (cfl * q / Np - Md)) * (M_PI * (cfl * q / Np - Md));
      return (1.0 - 2.0 * alpha) * std::exp(-alpha);
    }
  };

  template<typename T>
  struct ContinuousSource final : TemporalSource<T> {
    [[nodiscard]] T eval(const T t) const override {
      if (t < start_time or t > stop_time) { return 0.0; }
      return std::sin(omega * t - phase);
    }

    // T freq;
    T omega;
    T start_time;
    T stop_time;
    T phase;
    // T delay;
  };


  // // todo: Does this serve any purpose? Does it solve the soft-source amplitude problem? It's certainly overkill for a
  // //       simple spatial source to have an associated 1D solver...
  // //       Maybe replace with a lookup table?
  // template<typename T>
  // struct AuxiliarySource {
  //   using HIntegrator = FieldIntegrator1D<tf::types::Array1D<T>, FieldUpdate<Derivative::DX, Derivative::NoOp, true, size_t>>;
  //   using EIntegrator = FieldIntegrator1D<tf::types::Array1D<T>, FieldUpdate<Derivative::DX, Derivative::NoOp, false, size_t>>;
  //   using empty_t = tf::types::EmptyArray1D<T>;
  //
  //   static constexpr empty_t empty{};
  //
  //   static void updateH() { HIntegrator::apply(Hinc, Einc, empty, empty, Chyh, Chye, empty, empty, {0, 0, 0, 0, 0, 0}); }
  //   static void updateE() { EIntegrator::apply(Einc, Hinc, empty, empty, Ceze, Cezh, empty, empty, {1, 1, 1, 1, 0, 0}); }
  //   void apply_src(const T q) {
  //     auto result = 1.0;
  //     for (const auto& src : t_srcs) {
  //       result *= src->eval(q);
  //     }
  //     Einc[0] = result;
  //   }
  //
  //   std::vector<std::unique_ptr<TemporalSource<T>>> t_srcs{};
  //
  //   tf::types::Array1D<T> Einc;
  //   tf::types::Array1D<T> Ceze;
  //   tf::types::Array1D<T> Cezh;
  //
  //   tf::types::Array1D<T> Hinc;
  //   tf::types::Array1D<T> Chyh;
  //   tf::types::Array1D<T> Chye;
  // };



  template<typename T>
  struct SpatialSource {
    // SpatialSource(const size_t x0_, const size_t x1_)
    // requires (std::same_as<Array, tf::types::Array1D<T>>)
    // : x0(x0_), x1(x1_), y0(0), y1(1), z0(0), z1(1)
    // {}
    //
    // SpatialSource(const size_t x0_, const size_t x1_, const size_t y0_, const size_t y1_)
    // requires (std::same_as<Array, tf::types::Array2D<T>>)
    // : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(0), z1(1)
    // {}

    SpatialSource(const size_t x0_, const size_t x1_, const size_t y0_=0, const size_t y1_=1, const size_t z0_=0, const size_t z1_=1)
    : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_)
    {}

    [[nodiscard]] T eval(const T t) {
      auto result = 1.0;
      for (const auto& src : t_srcs) {
        result *= src->eval(t);
      }

      return result;
    }

    size_t x0, x1, y0, y1, z0, z1;
    std::vector<std::unique_ptr<TemporalSource<T>>> t_srcs{};
  };

  // todo: I can probably fold this into the SpatialSource class to get rid of a level of abstraction.
  template<typename Array>
  struct CurrentSource {
    using array_t = Array;
    using value_t = typename Array::value_t;
    using dimension_t = typename Array::dimension_t;

    CurrentSource(Array* f, SpatialSource<value_t>&& s)
    : field(f),
      src(std::move(s))
    {}

    // todo: These are soft sources, but don't fix the soft source problem (incorrect amplitude)
    //       Not sure what the best fix is. Auxilliary sources seem overkill. Maybe a lookup table?
    void apply(const value_t q) const
    requires (dimension_t::value == 1)
    {
      for (size_t i = src.x0; i <= src.x1; ++i) {
        (*field)[i] += src.eval(q);
      }
    }

    void apply(const value_t q)
    requires (dimension_t::value == 2)
    {
      for (size_t i = src.x0; i < src.x1; ++i) {
        for (size_t j = src.y0; j < src.y1; ++j) {
          (*field)(i, j) += src.eval(q);
        }
      }
    }

    void apply(const value_t q) const
    requires (dimension_t::value == 3)
    {
      for (size_t i = src.x0; i < src.x1; ++i) {
        for (size_t j = src.y0; j < src.y1; ++j) {
          for (size_t k = src.z0; k < src.z1; ++k) {
            (*field)(i, j, k) += src.eval(q);
          }
        }
      }
    }

    Array* const field;
    SpatialSource<value_t> src;
  };
  
  //
} // end namespace tf::electromagnetics::sources
#endif //EM_SOURCES_H
