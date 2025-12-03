#ifndef EM_SOURCES_HPP
#define EM_SOURCES_HPP

#include "em_params.hpp"
#include "program_params.hpp"
#include "constants.hpp"
#include "math_utils.hpp"
#include "mdspan.hpp"
#include "vec3.hpp"

#include <cassert>
#include <memory>
#include <print>
#include <utility>

namespace tf::electromagnetics {
struct TemporalSource {
   virtual ~TemporalSource() = default;
   [[nodiscard]] virtual double eval(double) const = 0;
};

struct RickerSource final : TemporalSource {
   explicit RickerSource(const double freq_)
   : freq(freq_) {}

   [[nodiscard]] double eval(const double t) const override {
      constexpr auto Md = 2.0;
      const auto alpha = math::SQR(constants::pi<double> * freq * (t - Md / freq));
      const auto temp = (1.0 - 2.0 * alpha) * std::exp(-alpha);
      return temp;
   }

   double freq;
}; // end struct RickerSource

struct GaussianSource final : TemporalSource {
   GaussianSource(const double width_, const double power_, const double delay_)
   : width(width_),
     power(power_),
     delay(delay_)
   {}

   [[nodiscard]] double eval(const double t) const override {
      constexpr auto tol = 1e-15;
      const auto     val = std::exp(-0.5 * std::pow((t - delay) / width, power));
      return val <= tol ? 0.0 : val;
   }

   double width;
   double power;
   double delay;
}; // end struct GaussianSource


struct ContinuousSource final : TemporalSource {
   explicit ContinuousSource(const double omega_, const double phase_, const double start_, const double stop_)
   : omega(omega_),
     start(start_),
     stop(stop_),
     phase(phase_)
   {}

   [[nodiscard]] double eval(const double t) const override {
      if (t < start or t > stop) { return 0.0; }
      return std::sin(omega * t - phase);
   }

   double omega;
   double start;
   double stop;
   double phase;
}; // end struct ContinuousSource

struct CurrentSource {
   using temporal_vec = std::vector<std::unique_ptr<TemporalSource>>;

   CurrentSource(const mdspan_t& f_, temporal_vec tsrc_, const double amp_)
   : field(f_),
     time_srcs(std::move(tsrc_)),
     amp(amp_)
   {}

   void apply(const double t) const {
      auto result = amp;
      for (const auto& s : time_srcs) { result *= s->eval(t); }
      for (auto i = 0zu; i < field.extent(0); i++) {
         for (auto j = 0zu; j < field.extent(1); j++) {
            for (auto k = 0zu; k < field.extent(2); k++) {
               field[i, j, k] = result;
            } // end for(k)
         } // end for(j)
      } // end for(i)
   } // end apply()

   mdspan_t field;
   temporal_vec time_srcs;
   double amp;
}; // end struct CurrentSource

struct GaussianBeam {
   using temporal_vec = std::vector<std::unique_ptr<TemporalSource>>;

   GaussianBeam(const mdspan_t& f_, temporal_vec tsrc_, const double amp_, const double omega_, const double w0_, const vec3<double>& wpos_)
   : field(f_),
     time_srcs(std::move(tsrc_)),
     coeffs(Nz - 2 * (BCDepth + 10)),
     waist_pos(wpos_),
     waist_size(w0_),
     amp(amp_)
   {
      constexpr auto x0 = static_cast<double>(BCDepth + 10);
      constexpr auto z0 = BCDepth + 10;
      constexpr auto z1 = Nz - z0;

      constexpr auto zmin = z_range[0] + dz * static_cast<double>(z0);
      constexpr auto zmax = z_range[0] + dz * static_cast<double>(z1 - 1);
      constexpr auto xpos = x_range[0] + (x0 * dx);

      const auto    z = waist_pos.x - xpos; // +x? direction
      const auto    k = omega_ / constants::c<double>;
      const auto   zR = 0.5 * k * math::SQR(waist_size);
      const auto   wz = waist_size * std::sqrt(1.0 + math::SQR(z / zR));
      const auto   RC = z * (1.0 + math::SQR(zR / z));
      const auto gouy = std::atan2(z, zR);
      const auto   c1 = waist_size / wz;

      const auto r = math::linspace(zmin, zmax, z1 - z0, true);
      const auto wz2 = wz * wz;
      for (std::size_t i = 0; i < r.size(); ++i) {
         const auto r2 = r[i] * r[i];
         coeffs[i] = c1 * std::exp(-r2 / wz2) * std::cos(0.5 * k * r2 / RC - gouy);
      }
   } // end GaussianBeam ctor

   void apply(const double t) const {
      auto result = amp;
      for (const auto& s : time_srcs) { result *= s->eval(t); }
      // if (result < 1.0e-10) { return; }

      for (auto i = 0zu; i < field.extent(0); i++) {
         for (auto j = 0zu; j < field.extent(1); j++) {
            for (auto k = 0zu; k < field.extent(2); k++) {
               field[i, j, k] += result;
            } // end for(k)
         } // end for(j)
      } // end for(i)
   } // end apply()

   mdspan_t field;
   temporal_vec time_srcs;
   std::vector<double> coeffs;
   vec3<double> waist_pos;
   double waist_size;
   double amp;
}; // end struct GaussianBeam

void add_gaussianbeam(auto& emdata) {
   using temporal_vec = std::vector<std::unique_ptr<TemporalSource>>;

   constexpr auto freq = constants::c<double> / 8.0e-7; // Hz -> c / 800 nm
   constexpr auto omega = 2.0 * constants::pi<double> * freq;

   constexpr auto amp = 1.583 * 2.75e13; // V/m
   constexpr auto w0 = 2.5479e-6; // meters, waste size

   constexpr auto width = 1.2739827e-14; // ~12.74 fs
   constexpr auto delay = 3.0 * width;

   constexpr vec3 waist_pos{0.0, 0.0, 0.0};

   auto make_srcvec = [&]() -> temporal_vec {
      temporal_vec result{};
      result.push_back(std::make_unique<GaussianSource>(width, 2.0, delay));
      result.push_back(std::make_unique<ContinuousSource>(omega, 0.0, 0.0, 1.0e30));
      return result;
   };

   emdata.beams.push_back(
      GaussianBeam{
      // mdspan_t{&emdata.Eyf[25, 0, 25], {std::extents{1, 1, Nz - 50}, ey_stride}},
      mdspan_t{&emdata.Eyf[BCDepth + 10, 0, BCDepth + 10], {std::extents{1, Ny - 1, Nz - 2 * (BCDepth + 10)}, ey_stride}},
      make_srcvec(),
      amp,
      omega,
      w0,
      waist_pos}
   );
} // end add_gaussianbeam()
} // end namespace tf::electromagnetics

#endif //EM_SOURCES_HPP
