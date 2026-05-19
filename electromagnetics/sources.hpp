#ifndef EM_SOURCES_HPP
#define EM_SOURCES_HPP

#include "program_params.hpp"
#include "constants.hpp"
#include "array.hpp"
#include "math_utils.hpp"

// #include <cassert>
#include <memory>
#include <print>

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
      const auto alpha = math::SQR(static_cast<double>(constants::pi) * freq * (t - Md / freq));
      const auto temp = (1.0 - 2.0 * alpha) * std::exp(-alpha);
      return temp;
   }

   double freq;
}; // end struct RickerSource

struct BlackmanHarris final : TemporalSource {
   /*
    * https://en.wikipedia.org/wiki/Window_function#Blackman%E2%80%93Harris_window
    * Blackman-Harris Window is similar to a Gaussian except having better frequency
    * behaviors to avoid more DC components.
    *
    * cutoff_factor is equivalent to N/2 -> 1.55 / (2*freq) -> 1.55 / (omega/pi) -> (1.55*pi)/omega -> 4.867/omega
    * bn is from 2*pi*t/N -> 2*pi*t / (1.55 / freq) -> 2*pi*freq*t / 1.55 -> omega * t / 1.55,
    *    so the array is {1/1.55, 2/1.55, 3/1.55}
    *
    */
   explicit BlackmanHarris(const double dx_)
   : omega((constants::pi2 * constants::c) / (Nl * dx_)),
     duration{cutoff_factor / omega}
   {}

   [[nodiscard]] double eval(const double t) const override {
      if (t > duration) { return 1.0; }

      const auto c1 = std::cos(bn[0] * omega * t);
      const auto c2 = std::cos(bn[1] * omega * t);
      const auto c3 = std::cos(bn[2] * omega * t);
      return an[0] + (an[1] * c1) + (an[2] * c2) + (an[3] * c3);
   }

   double                  omega;
   double                  duration;
   static constexpr double Nl = 25.0;
   static constexpr double cutoff_factor{4.867f};
   static constexpr double an[4]{0.3588, -0.4882, 0.1413, -0.0116};
   static constexpr double bn[3]{0.6452, 1.2903, 1.9355};
}; // end struct BlackmanHarris


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

struct SpatialSource {
   using temporal_vec = std::vector<std::unique_ptr<TemporalSource>>;
   using offset_t     = std::array<std::size_t, 6>;

   SpatialSource(temporal_vec&& ts, const double amp_, const offset_t& o_)
   : amplitude(amp_),
     offsets{o_},
     t_srcs{std::move(ts)}
   {}

   [[nodiscard]] double eval(const double t) const {
      auto result = amplitude;
      for (const auto& src: t_srcs) {
         result *= src->eval(t);
      }
      return result;
   } // end SpatialSource::eval

   double amplitude;
   offset_t offsets;
   temporal_vec t_srcs;
}; // end struct SpatialSource

struct CurrentSource {
   using array_t = Array3D<double>;

   CurrentSource(array_t* const f, SpatialSource&& s)
   : field(f),
     src(std::move(s))
   {}

   void apply(const double t) const {
      const auto& [x0, x1, y0, y1, z0, z1] = src.offsets;
      const auto  val                      = src.eval(t);
      if (val == 0.0) { return; }

      for (size_t i = x0; i < x1; ++i) {
         for (size_t j = y0; j < y1; ++j) {
            for (size_t k = z0; k < z1; ++k) {
               (*field)(i, j, k) = val;
            }
         }
      }
   }

   array_t* const field;
   SpatialSource  src;
}; // end struct CurrentSource

struct GaussianBeam {
   using array_t  = Array3D<double>;
   using offset_t = std::array<std::size_t, 6>;

   explicit GaussianBeam(array_t& field_)
   : field(field_),
     zs(math::linspace(z_range[0], z_range[1], Nz - 1))
   {}

   void apply(const double t) const {
      if (t > 60.0e-15) { return; }

      constexpr auto x0 = PMLDepth + 10zu;
      constexpr auto z0 = PMLDepth + 10zu;
      constexpr auto z1 = Nz - z0 - 1zu;

      constexpr auto lambda = 8.0e-7;
      constexpr auto omega = 2.0 * constants::pi * constants::c / lambda;
      constexpr auto omega_env = constants::pi / 60.0e-15;
      constexpr auto E0 = -2.75e13; // V/m
      constexpr auto w0 = 2.5479e-6;   // meters, waste size
      constexpr auto xspot = 15.0e-6;
      constexpr auto xR = constants::pi * math::SQR(w0) / lambda;
      constexpr auto RC = xspot * (1.0 + math::SQR(xR / xspot));
      constexpr auto kn = 2.0 * constants::pi / lambda;

      const auto wx = w0 * std::sqrt(1.0 + math::SQR(xspot / xR));
      const auto gouy = std::atan(xspot / xR);
      const auto c1 = 1.288 * E0 * w0 / wx; // Fudge it, fudge it all

      for (auto k = z0; k < z1; ++k) {
         const auto kdx = k - z0;
         field(x0, 0, k) += c1 * std::exp(-math::SQR(zs[kdx] / wx))
                               * std::sin(omega_env * t)
                               * std::sin(omega * t + 0.5 * kn * math::SQR(zs[kdx]) / RC - gouy);
      }
   }

   array_t& field;
   std::vector<double> zs;
};

} // end namespace tf::electromagnetics

#endif //EM_SOURCES_HPP
