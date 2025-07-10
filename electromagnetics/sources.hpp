#ifndef EM_SOURCES_HPP
#define EM_SOURCES_HPP

#include "program_params.hpp"
#include "constants.hpp"
#include "array.hpp"
#include "math_utils.hpp"

#include "dbg.h"

#include <memory>
#include <vector>
#include <cassert>
// #include <fstream>
// #include <iomanip>

namespace tf::electromagnetics {
struct TemporalSource {
   virtual ~TemporalSource() = default;
   [[nodiscard]] virtual compute_t eval(compute_t) const = 0;
};


struct RickerSource final : TemporalSource {
   explicit RickerSource(const compute_t freq_)
   : freq(freq_) {}

   [[nodiscard]] compute_t eval(const compute_t t) const override {
      constexpr auto Md = 2.0_fp;
      const auto alpha = math::SQR(static_cast<compute_t>(constants::pi<compute_t>) * freq * (t - Md / freq));
      return (1.0_fp - 2.0_fp * alpha) * std::exp(-alpha);
   }

   compute_t freq;
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
   explicit BlackmanHarris(const compute_t dx_)
   : omega((constants::pi2<compute_t> * constants::c<compute_t>) / (Nl * dx_)),
     duration{cutoff_factor / omega}
   {}

   [[nodiscard]] compute_t eval(const compute_t t) const override {
      if (t > duration) { return 1.0_fp; }

      const auto c1 = std::cos(bn[0] * omega * t);
      const auto c2 = std::cos(bn[1] * omega * t);
      const auto c3 = std::cos(bn[2] * omega * t);
      return an[0] + (an[1] * c1) + (an[2] * c2) + (an[3] * c3);
   }

   compute_t                  omega;
   compute_t                  duration;
   static constexpr compute_t Nl = 25.0_fp;
   static constexpr compute_t cutoff_factor{4.867f};
   static constexpr compute_t an[4]{0.3588_fp, -0.4882_fp, 0.1413_fp, -0.0116_fp};
   static constexpr compute_t bn[3]{0.6452_fp, 1.2903_fp, 1.9355_fp};
}; // end struct BlackmanHarris


struct GaussianSource final : TemporalSource {
   GaussianSource(const compute_t width_, const compute_t power_, const compute_t delay_)
   : width(width_),
     power(power_),
     delay(delay_)
   {}

   [[nodiscard]] compute_t eval(const compute_t t) const override {
      constexpr auto tol = 1e-15_fp;
      const auto     val = std::exp(-0.5_fp * std::pow((t - delay) / width, power));
      return val <= tol ? 0.0_fp : val;
   }

   compute_t width;
   compute_t power;
   compute_t delay;
}; // end struct GaussianSource


struct ContinuousSource final : TemporalSource {
   explicit ContinuousSource(const compute_t omega_, const compute_t phase_, const compute_t start_,
                             const compute_t stop_, const compute_t)
   : omega(omega_),
     start(start_),
     stop(stop_),
     phase(phase_) //, ramp{dx_}
   {}

   [[nodiscard]] compute_t eval(const compute_t t) const override {
      if (t < start or t > stop) { return 0.0_fp; }
      return std::sin(omega * t - phase);
   }

   compute_t omega;
   compute_t start;
   compute_t stop;
   compute_t phase;
   // BlackmanHarris ramp;
}; // end struct ContinuousSource

struct SpatialSource {
   using temporal_vec = std::vector<std::unique_ptr<TemporalSource>>;
   using offset_t     = std::array<std::size_t, 6>;

   SpatialSource(temporal_vec&& ts, const compute_t amp_, const offset_t& o_)
   : amplitude(amp_),
     offsets{o_},
     t_srcs{std::move(ts)}
   {}

   [[nodiscard]] compute_t eval(const compute_t t) const {
      auto result = amplitude;
      for (const auto& src: t_srcs) {
         result *= src->eval(t);
      }
      return result;
   } // end SpatialSource::eval

   compute_t    amplitude;
   offset_t     offsets;
   temporal_vec t_srcs;
}; // end struct SpatialSource

struct CurrentSource {
   using array_t = Array3D<compute_t>;

   CurrentSource(array_t* const f, SpatialSource&& s)
   : field(f),
     src(std::move(s))
   {}

   // void apply(compute_t) const;

   array_t* const field;
   SpatialSource  src;
}; // end struct CurrentSource

struct GaussianBeam : CurrentSource {
   using array_t  = Array3D<compute_t>;
   using offset_t = std::array<std::size_t, 6>;

   GaussianBeam(array_t* const         f_,
                const compute_t        waist_,
                const compute_t        omega_,
                const vec3<compute_t>& waist_pos_,
                SpatialSource&&        s_)
   : CurrentSource(f_, std::forward<SpatialSource>(s_)),
     waist_size(waist_),
     waist_pos(waist_pos_),
     coeffs(src.offsets[5] - src.offsets[4])
   {
      const auto& [x0, x1, y0, y1, z0, z1] = src.offsets;
      assert((z1 - z0) == coeffs.size());
      const auto xpos = x_range[0] + (static_cast<compute_t>(x0) * dx);
      const auto z    = xpos - waist_pos[0]; // -x direction
      // const auto z    = waist_pos[0] - xpos; // +x? direction

      assert(z != 0.0_fp);
      const auto k  = omega_ / static_cast<compute_t>(constants::c<compute_t>);
      const auto zR = 0.5_fp * k * math::SQR(waist_size);
      const auto wz = waist_size * std::sqrt(1.0_fp + math::SQR(z / zR));
      const auto RC   = z * (1.0_fp + math::SQR(zR / z));
      const auto gouy = std::atan2(z, zR);
      const auto c1   = waist_size / wz;
      const auto zmin = z_range[0] + dz * static_cast<compute_t>(z0);
      const auto zmax = z_range[0] + dz * static_cast<compute_t>(z1 - 1);
      const auto r = math::linspace(zmin, zmax, z1 - z0, true);
      const auto wz2 = wz * wz;
      for (std::size_t i = 0; i < r.size(); ++i) {
         const auto r2 = r[i] * r[i];
         coeffs[i] = c1 * std::exp(-r2 / wz2) * std::cos(0.5_fp * k * r2 / RC - gouy);
      }
   } // end GaussianBeam ctor

   void apply(const compute_t t) const {
      const auto& [x0, x1, y0, y1, z0, z1] = src.offsets;
      const auto  val                      = src.eval(t);
      if (val == 0.0_fp) { return; }

      for (size_t i = x0; i < x1; ++i) {
         for (size_t j = y0; j < y1; ++j) {
            for (size_t k = z0; k < z1; ++k) {
               (*field)(i, j, k) += coeffs[k - z0] * val;
            }
         }
      }
   }

   compute_t              waist_size;
   vec3<compute_t>        waist_pos;
   std::vector<compute_t> coeffs;
}; // end struct GaussianBeam

void add_gaussianbeam(auto& em) {
   using temporal_vec = std::vector<std::unique_ptr<TemporalSource>>;

   constexpr auto freq = constants::c<compute_t> / 8.0e-7_fp; // Hz -> c / 800 nm
   constexpr auto omega = 2.0_fp * constants::pi<compute_t> * freq;

   constexpr auto amp = 1.1750978 * 2.75e13_fp; // V/m
   constexpr auto w0 = 2.5479e-6_fp; // meters, waste size

   constexpr auto width = 1.2739827e-14_fp; // ~12.74 fs
   constexpr auto delay = 2.0 * width;

   vec3 waist_pos{0.0_fp, 0.0_fp, 0.0_fp};

   constexpr auto x0 = PMLDepth + 20zu;
   constexpr auto x1 = x0 + 1;
   constexpr auto y0 = 0zu;
   constexpr auto y1 = 1zu;
   constexpr auto z0 = PMLDepth + 20zu;
   constexpr auto z1 = Nz - z0;

   using continuous_t = ContinuousSource;
   auto make_continuous = [&](temporal_vec& srcs) {
      srcs.push_back(std::make_unique<continuous_t>(omega, 0.0f, 0.0f, 1.0e30f, dx));
   };

   using gaussian_t = GaussianSource;
   auto make_gaussian = [&](temporal_vec& srcs) {
      srcs.push_back(std::make_unique<gaussian_t>(width, 2.0_fp, delay));
   };


   auto make_srcvec = [&]() -> temporal_vec {
      temporal_vec result{};
      make_gaussian(result);
      make_continuous(result);
      return result;
   };

   em.emdata.srcs.emplace_back(
      &em.emdata.Ey,
      w0,
      omega,
      waist_pos,
      SpatialSource(
         make_srcvec(),
         amp,
         {x0, x1, y0, y1, z0, z1}
      )
   );
}
} // end namespace tf::electromagnetics

#endif //EM_SOURCES_HPP
