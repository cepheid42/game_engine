#ifndef EM_SOURCES_HPP
#define EM_SOURCES_HPP

#include "program_params.hpp"
#include "constants.hpp"
#include "array.hpp"
#include "math_utils.hpp"

#include <memory>

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

struct GaussianBeam : CurrentSource {
   using array_t  = Array3D<double>;
   using offset_t = std::array<std::size_t, 6>;

   GaussianBeam(array_t* const         f_,
                const double        waist_,
                const double        omega_,
                const vec3<double>& waist_pos_,
                SpatialSource&&        s_)
   : CurrentSource(f_, std::forward<SpatialSource>(s_)),
     waist_size(waist_),
     waist_pos(waist_pos_),
     coeffs(src.offsets[5] - src.offsets[4])
   {
      const auto& [x0, x1, y0, y1, z0, z1] = src.offsets;
      assert((z1 - z0) == coeffs.size());
      const auto xpos = x_range[0] + (static_cast<double>(x0) * dx);
      const auto z    = waist_pos[0] - xpos; // +x? direction
      assert(z != 0.0);
      const auto k  = omega_ / constants::c;
      const auto zR = 0.5 * k * math::SQR(waist_size);
      const auto wz = waist_size * std::sqrt(1.0 + math::SQR(z / zR));
      const auto RC   = z * (1.0 + math::SQR(zR / z));
      const auto gouy = std::atan2(z, zR);
      const auto c1   = waist_size / wz;
      const auto zmin = z_range[0] + dz * static_cast<double>(z0);
      const auto zmax = z_range[0] + dz * static_cast<double>(z1 - 1);
      const auto r = math::linspace(zmin, zmax, z1 - z0, true);
      const auto wz2 = wz * wz;
      for (std::size_t i = 0; i < r.size(); ++i) {
         const auto r2 = r[i] * r[i];
         coeffs[i] = c1 * std::exp(-r2 / wz2) * std::cos(0.5 * k * r2 / RC - gouy);
      }
   } // end GaussianBeam ctor

   void apply(const double t) const {
      const auto& [x0, x1, y0, y1, z0, z1] = src.offsets;
      const auto  val                      = src.eval(t);
      if (val == 0.0) { return; }

      for (size_t i = x0; i < x1; ++i) {
         for (size_t j = y0; j < y1; ++j) {
            for (size_t k = z0; k < z1; ++k) {
               (*field)(i, j, k) += coeffs[k - z0] * val;
            }
         }
      }
   }

   double              waist_size;
   vec3<double>        waist_pos;
   std::vector<double> coeffs;
}; // end struct GaussianBeam

void add_linesource(auto& emsolver, const std::string component, const double amp, const std::array<double,6> bounds, const double omega, const double phase = 0.0, const double start = 0.0, const double stop = 1.0e30) {

   using temporal_vec = std::vector<std::unique_ptr<TemporalSource>>;
   using continuous_t = ContinuousSource;

   // convert real-space bounds to node offsets
   std::array offsets = {
      static_cast<std::size_t>(std::floor((bounds[0] - x_range[0]) / dx)),
      static_cast<std::size_t>(std::floor((bounds[1] - x_range[0]) / dx)),
      static_cast<std::size_t>(std::floor((bounds[2] - y_range[0]) / dy)),
      static_cast<std::size_t>(std::floor((bounds[3] - y_range[0]) / dy)),
      static_cast<std::size_t>(std::floor((bounds[4] - z_range[0]) / dz)),
      static_cast<std::size_t>(std::floor((bounds[5] - z_range[0]) / dz))
   };

   // ensure the line of nodes is at least one node thick
   if (offsets[0] == offsets[1]) { offsets[1]++; }
   if (offsets[2] == offsets[3]) { offsets[3]++; }
   if (offsets[4] == offsets[5]) { offsets[5]++; }

   // create temporal source
   temporal_vec result{};
   result.push_back(std::make_unique<continuous_t>(omega, phase, start, stop));

   // create spatial source
   if (component == "ex") {
      emsolver.emdata.srcs.emplace_back(&emsolver.emdata.Ex, SpatialSource(std::move(result), amp, offsets));
   }
   else if (component == "ey") {
      emsolver.emdata.srcs.emplace_back(&emsolver.emdata.Ey, SpatialSource(std::move(result), amp, offsets));
   }
   else if (component == "ez") {
      emsolver.emdata.srcs.emplace_back(&emsolver.emdata.Ez, SpatialSource(std::move(result), amp, offsets));
   }
   else {
      throw std::invalid_argument("add_linesource: Field Component can only be ['ex', 'ey', 'ez'] got '" + component + "' instead.");
   }

} // end add_linesource


void add_rmf_antennas(auto& emsolver, const auto rmf_params) {
   const auto [amp, freq, antenna_lz, antenna_lxy, antenna_r] = rmf_params;
   const auto period = 1.0 / freq;
   const auto omega = 2 * constants::pi * freq;
   const auto delay_lr = period / 4.0;

   const auto amp_m = -amp;
   const auto amp2 = 2 * amp;

   // Segments are numbered 1-7 counterclockwise, starting at z0 side. Segment 7 is the centerline
   // The "left" antenna is the z0 side, "right" is z1 side.
   //    6   5
   //   --- ---
   // 1|  7|   | 4
   //   --- ---
   //    2   3

   // antenna origin is in center of center conductor
   const std::array origin_t = {0.0, antenna_r, 0.0};
   const std::array origin_b = {0.0, -antenna_r, 0.0};
   const std::array origin_l = {-antenna_r, 0.0, 0.0};
   const std::array origin_r = {antenna_r, 0.0, 0.0};

   const auto phase_tb = -constants::pi / 2;
   const auto start_tb = 0.0;
   const auto phase_lr = 0.0;
   const auto start_lr = delay_lr;

   // bottom (-y) antenna
   add_linesource(emsolver, "ex", amp_m, {origin_b[0] - antenna_lxy / 2.0, origin_b[0] + antenna_lxy / 2.0, origin_b[1], origin_b[1], origin_b[2] - antenna_lz / 2.0, origin_b[2] - antenna_lz / 2.0}, omega, phase_tb, start_tb);
   add_linesource(emsolver, "ex", amp_m, {origin_b[0] - antenna_lxy / 2.0, origin_b[0] + antenna_lxy / 2.0, origin_b[1], origin_b[1], origin_b[2] + antenna_lz / 2.0, origin_b[2] + antenna_lz / 2.0}, omega, phase_tb, start_tb);

   add_linesource(emsolver, "ez", amp, {origin_b[0] - antenna_lxy / 2.0, origin_b[0] - antenna_lxy / 2.0, origin_b[1], origin_b[1], origin_b[2] - antenna_lz / 2.0, origin_b[2]}, omega, phase_tb, start_tb);
   add_linesource(emsolver, "ez", amp_m, {origin_b[0] + antenna_lxy / 2.0, origin_b[0] + antenna_lxy / 2.0, origin_b[1], origin_b[1], origin_b[2] - antenna_lz / 2.0, origin_b[2]}, omega, phase_tb, start_tb);

   add_linesource(emsolver, "ez", amp_m, {origin_b[0] - antenna_lxy / 2.0, origin_b[0] - antenna_lxy / 2.0, origin_b[1], origin_b[1], origin_b[2], origin_b[2] + antenna_lz / 2.0}, omega, phase_tb, start_tb);
   add_linesource(emsolver, "ez", amp, {origin_b[0] + antenna_lxy / 2.0, origin_b[0] + antenna_lxy / 2.0, origin_b[1], origin_b[1], origin_b[2], origin_b[2] + antenna_lz / 2.0}, omega, phase_tb, start_tb);

   add_linesource(emsolver, "ex", amp2, {origin_b[0] - antenna_lxy / 2.0, origin_b[0] + antenna_lxy / 2.0, origin_b[1], origin_b[1], origin_b[2], origin_b[2]}, omega, phase_tb, start_tb);

   // top (+y) antenna
   add_linesource(emsolver, "ex", amp_m, {origin_t[0] - antenna_lxy / 2.0, origin_t[0] + antenna_lxy / 2.0, origin_t[1], origin_t[1], origin_t[2] - antenna_lz / 2.0, origin_t[2] - antenna_lz / 2.0}, omega, phase_tb, start_tb);
   add_linesource(emsolver, "ex", amp_m, {origin_t[0] - antenna_lxy / 2.0, origin_t[0] + antenna_lxy / 2.0, origin_t[1], origin_t[1], origin_t[2] + antenna_lz / 2.0, origin_t[2] + antenna_lz / 2.0}, omega, phase_tb, start_tb);

   add_linesource(emsolver, "ez", amp, {origin_t[0] - antenna_lxy / 2.0, origin_t[0] - antenna_lxy / 2.0, origin_t[1], origin_t[1], origin_t[2] - antenna_lz / 2.0, origin_t[2]}, omega, phase_tb, start_tb);
   add_linesource(emsolver, "ez", amp_m, {origin_t[0] + antenna_lxy / 2.0, origin_t[0] + antenna_lxy / 2.0, origin_t[1], origin_t[1], origin_t[2] - antenna_lz / 2.0, origin_t[2]}, omega, phase_tb, start_tb);

   add_linesource(emsolver, "ez", amp_m, {origin_t[0] - antenna_lxy / 2.0, origin_t[0] - antenna_lxy / 2.0, origin_t[1], origin_t[1], origin_t[2], origin_t[2] + antenna_lz / 2.0}, omega, phase_tb, start_tb);
   add_linesource(emsolver, "ez", amp, {origin_t[0] + antenna_lxy / 2.0, origin_t[0] + antenna_lxy / 2.0, origin_t[1], origin_t[1], origin_t[2], origin_t[2] + antenna_lz / 2.0}, omega, phase_tb, start_tb);

   add_linesource(emsolver, "ex", amp2, {origin_t[0] - antenna_lxy / 2.0, origin_t[0] + antenna_lxy / 2.0, origin_t[1], origin_t[1], origin_t[2], origin_t[2]}, omega, phase_tb, start_tb);

   // left (-x) antenna
   add_linesource(emsolver, "ey", amp_m, {origin_l[0], origin_l[0], origin_l[1] - antenna_lxy / 2.0, origin_l[1] + antenna_lxy / 2.0, origin_l[2] - antenna_lz / 2.0, origin_l[2] - antenna_lz / 2.0}, omega, phase_lr, start_lr);
   add_linesource(emsolver, "ey", amp_m, {origin_l[0], origin_l[0], origin_l[1] - antenna_lxy / 2.0, origin_l[1] + antenna_lxy / 2.0, origin_l[2] + antenna_lz / 2.0, origin_l[2] + antenna_lz / 2.0}, omega, phase_lr, start_lr);

   add_linesource(emsolver, "ez", amp, {origin_l[0], origin_l[0], origin_l[1] - antenna_lxy / 2.0, origin_l[1] - antenna_lxy / 2.0, origin_l[2] - antenna_lz / 2.0, origin_l[2]}, omega, phase_lr, start_lr);
   add_linesource(emsolver, "ez", amp_m, {origin_l[0], origin_l[0], origin_l[1] + antenna_lxy / 2.0, origin_l[1] + antenna_lxy / 2.0, origin_l[2] - antenna_lz / 2.0, origin_l[2]}, omega, phase_lr, start_lr);

   add_linesource(emsolver, "ez", amp_m, {origin_l[0], origin_l[0], origin_l[1] - antenna_lxy / 2.0, origin_l[1] - antenna_lxy / 2.0, origin_l[2], origin_l[2] + antenna_lz / 2.0}, omega, phase_lr, start_lr);
   add_linesource(emsolver, "ez", amp, {origin_l[0], origin_l[0] , origin_l[1] + antenna_lxy / 2.0, origin_l[1]+ antenna_lxy / 2.0, origin_l[2], origin_l[2] + antenna_lz / 2.0}, omega, phase_lr, start_lr);

   add_linesource(emsolver, "ey", amp2, {origin_l[0], origin_l[0], origin_l[1] - antenna_lxy / 2.0, origin_l[1] + antenna_lxy / 2.0, origin_l[2], origin_l[2]}, omega, phase_lr, start_lr);

   // right (+x) antenna
   add_linesource(emsolver, "ey", amp_m, {origin_r[0], origin_r[0], origin_r[1] - antenna_lxy / 2.0, origin_r[1] + antenna_lxy / 2.0, origin_r[2] - antenna_lz / 2.0, origin_r[2] - antenna_lz / 2.0}, omega, phase_lr, start_lr);
   add_linesource(emsolver, "ey", amp_m, {origin_r[0], origin_r[0], origin_r[1] - antenna_lxy / 2.0, origin_r[1] + antenna_lxy / 2.0, origin_r[2] + antenna_lz / 2.0, origin_r[2] + antenna_lz / 2.0}, omega, phase_lr, start_lr);

   add_linesource(emsolver, "ez", amp, {origin_r[0], origin_r[0], origin_r[1] - antenna_lxy / 2.0, origin_r[1] - antenna_lxy / 2.0, origin_r[2] - antenna_lz / 2.0, origin_r[2]}, omega, phase_lr, start_lr);
   add_linesource(emsolver, "ez", amp_m, {origin_r[0], origin_r[0], origin_r[1] + antenna_lxy / 2.0, origin_r[1] + antenna_lxy / 2.0, origin_r[2] - antenna_lz / 2.0, origin_r[2]}, omega, phase_lr, start_lr);

   add_linesource(emsolver, "ez", amp_m, {origin_r[0], origin_r[0], origin_r[1] - antenna_lxy / 2.0, origin_r[1] - antenna_lxy / 2.0, origin_r[2], origin_r[2] + antenna_lz / 2.0}, omega, phase_lr, start_lr);
   add_linesource(emsolver, "ez", amp, {origin_r[0], origin_r[0] , origin_r[1] + antenna_lxy / 2.0, origin_r[1]+ antenna_lxy / 2.0, origin_r[2], origin_r[2] + antenna_lz / 2.0}, omega, phase_lr, start_lr);

   add_linesource(emsolver, "ey", amp2, {origin_r[0], origin_r[0], origin_r[1] - antenna_lxy / 2.0, origin_r[1] + antenna_lxy / 2.0, origin_r[2], origin_r[2]}, omega, phase_lr, start_lr);
}

void add_gaussianbeam(auto& em) {
   using temporal_vec = std::vector<std::unique_ptr<TemporalSource>>;

   constexpr auto freq = constants::c / 8.0e-7; // Hz -> c / 800 nm
   constexpr auto omega = 2.0 * constants::pi * freq;

   constexpr auto amp = 1.583 * 2.75e13; // V/m
   constexpr auto w0 = 2.5479e-6; // meters, waste size

   constexpr auto width = 1.2739827e-14; // ~12.74 fs
   constexpr auto delay = 3.0 * width;

   vec3 waist_pos{0.0, 0.0, 0.0};

   constexpr auto x0 = PMLDepth + 10lu;
   constexpr auto x1 = x0 + 1;
   constexpr auto y0 = 0lu;
   constexpr auto y1 = 1lu;
   constexpr auto z0 = PMLDepth + 10lu;
   constexpr auto z1 = Nz - z0;

   using continuous_t = ContinuousSource;
   auto make_continuous = [&](temporal_vec& srcs) {
      srcs.push_back(std::make_unique<continuous_t>(omega, 0.0, 0.0, 1.0e30));
   };

   using gaussian_t = GaussianSource;
   auto make_gaussian = [&](temporal_vec& srcs) {
      srcs.push_back(std::make_unique<gaussian_t>(width, 2.0, delay));
   };

   auto make_srcvec = [&]() -> temporal_vec {
      temporal_vec result{};
      make_gaussian(result);
      make_continuous(result);
      return result;
   };

   em.emdata.beams.emplace_back(
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
