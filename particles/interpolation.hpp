#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

#include "math_utils.hpp"
#include "program_params.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "constants.hpp"

namespace tf::interp
{
/*---------------------------------------------------------------/
/-                     Interpolation Classes                    -/
/---------------------------------------------------------------*/
struct NGP {
   static constexpr std::size_t Begin   = 0;
   static constexpr std::size_t End     = 0;
   static constexpr std::size_t Order   = 0;
   static constexpr std::size_t Support = 1;

   static constexpr auto eval(const double) { return 1.0; }
   static constexpr auto shape_array(const double) { return std::array{1.0}; }
   static constexpr auto ds_array(const auto, const auto&) { return std::array{1.0}; }
};

struct CIC {
   static constexpr std::size_t Begin   = 0;
   static constexpr std::size_t End     = 1;
   static constexpr std::size_t Order   = 1;
   static constexpr std::size_t Support = 2;

   static constexpr auto eval(const auto x) { return 1.0 - std::abs(x); }

   static constexpr auto shape_array(const auto x) {
      return std::array{
         eval(x),
         eval(x - 1.0)
      };
   }

   static constexpr auto ds_array(const auto x1, const auto& s0) {
      return std::array{
         eval(x1)       - s0[0],
         eval(x1 - 1.0) - s0[1]
      };
   }
};

struct TSC {
   static constexpr std::size_t Begin   = 0;
   static constexpr std::size_t End     = 2;
   static constexpr std::size_t Order   = 2;
   static constexpr std::size_t Support = 3;

   static constexpr auto innerRadius(const auto x) { return 0.75 - math::SQR(x); }
   static constexpr auto outerRadius(const auto x) { return 0.5 * math::SQR(1.5 - x); }

   static constexpr auto eval(const auto x) {
      const auto absx = std::abs(x);
      return absx <= 0.5 ? innerRadius(absx) : outerRadius(absx);
   }

   static constexpr auto shape_array(const auto x) {
      return std::array{
         eval(x),
         eval(x - 1.0), // todo: Are these shifted correctly?
         eval(x - 2.0)
      };
   }

   static constexpr auto ds_array(const auto x1, const auto& s0) {
      return std::array{
         eval(x1)       - s0[0],
         eval(x1 - 1.0) - s0[1], // todo: Are these shifted correctly?
         eval(x1 - 2.0) - s0[2]
      };
   }
};

template<int ShapeOrder> struct InterpolationShape;
template<> struct InterpolationShape<0> { using Type = NGP; };
template<> struct InterpolationShape<1> { using Type = CIC; };
template<> struct InterpolationShape<2> { using Type = TSC; };

template<int D, typename Outer, typename Middle, typename Inner>
struct InterpolationStrategy {
   using OuterShape = Outer;
   using MiddleShape = Middle;
   using InnerShape = Inner;
   static constexpr auto permute(auto, auto, auto) = delete;
};

template<typename Outer, typename Middle, typename Inner>
struct InterpolationStrategy<0, Outer, Middle, Inner> {
   using OuterShape = Outer;
   using MiddleShape = Middle;
   using InnerShape = Inner;
   static constexpr auto permute(const auto x, const auto y, const auto z) { return std::tuple{z, x, y}; }
};

template<typename Outer, typename Middle, typename Inner>
struct InterpolationStrategy<1, Outer, Middle, Inner> {
   using OuterShape = Outer;
   using MiddleShape = Middle;
   using InnerShape = Inner;
   static constexpr auto permute(const auto x, const auto y, const auto z) { return std::tuple{y, z, x}; }
};

template<typename Outer, typename Middle, typename Inner>
struct InterpolationStrategy<2, Outer, Middle, Inner> {
   using OuterShape = Outer;
   using MiddleShape = Middle;
   using InnerShape = Inner;
   static constexpr auto permute(const auto x, const auto y, const auto z) { return std::tuple{x, y, z}; }
};


/*---------------------------------------------------------------/
/-                     Table Interpolation                      -/
/---------------------------------------------------------------*/
struct Table {
   std::vector<double> xs{};
   std::vector<double> ys{};

   explicit Table() = default;

   explicit Table(const std::string& filepath, const auto xscale=1.0, const auto yscale=1.0) {
      std::ifstream file(filepath);

      if (!file.is_open()) {
         throw std::runtime_error("Table: Unable to open cross section table: " + filepath);
      }

      double x{};
      double y{};

      std::string line;
      // std::getline(file, line);
      std::istringstream buff(line);
      while (std::getline(file, line)) {
         std::istringstream buffer(line);
         buffer >> x >> y;
         xs.push_back(x * xscale);
         ys.push_back(y * yscale);
      }
      file.close();
   }

   bool is_outofbounds(const auto e) const {
      return e <= xs[0] or e > xs[xs.size() - 1];
   }

   auto lerp(const auto e) const {
      const auto len = xs.size() - 1;
      // todo: not convinced these checks do anything useful
      if (e < xs[1]) { return ys[0]; }
      if (e > xs[len - 1]) { return ys[len]; }
      const auto upper = std::ranges::upper_bound(xs, e);
      const auto idx = std::ranges::distance(xs.cbegin(), upper) - 1;
      const auto slope = (e - xs[idx]) / (xs[idx + 1] - xs[idx]);
      return std::lerp(ys[idx], ys[idx + 1], slope);
   }
};

struct MultiTable {
   std::vector<std::vector<double>> data;
   std::array<double, 14> SB_k_over_gm1 = {1.0e-7, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.97, 0.99, 1.0};

   explicit MultiTable() = default;
   explicit MultiTable(const std::string& filepath)
   : data(16)
   {
      using namespace tf;
      std::ifstream file(filepath);

      if (!file.is_open()) {
         throw std::runtime_error("MultiTable: Unable to open cross section table: " + filepath);
      }

      double val{};
      std::string line;
      while (std::getline(file, line)) {
         std::istringstream buffer(line);
         auto col{0zu};
         while (buffer >> val) {
            data[col].push_back(val);
            col++;
         }
      }
      file.close();

      const auto n_energies = data[0].size();
      const auto n_columns = SB_k_over_gm1.size();

      for (auto i = 0zu; i < n_energies; ++i) {
         data[0][i] *= 1.0e6; // MeV -> eV

         // The input columns kE1...kE14 are k*d(sigma)/dk evaluated at the
         // tabulated values of x = k/(gamma-1).  The photon-number cross
         // section is therefore the integral of k*d(sigma)/dk over dln(k),
         // not a linear integral over x or k.
         //
         // The previous code created an artificial first interval 0 -> x_min,
         // shifted all subsequent intervals by one column, and then sampled
         // linearly in x.  That produces many photons below the nominal cutoff
         // and too few photons in the low-energy logarithmic bins.

         // Save the raw k*d(sigma)/dk data before overwriting the table columns
         // with the normalized CDF values.
         std::array<double, 14> kdsdk{};
         for (auto j = 0zu; j < n_columns; ++j) {
            kdsdk[j] = data[j + 1][i];
         }

         // Build the cumulative distribution beginning at x_min.  No photons
         // are sampled below SB_k_over_gm1[0].  The first CDF value is zero at
         // x_min, and the last CDF value is one at x_max.
         data[1][i] = 0.0;
         double sigma_ttl = 0.0;

         for (auto j = 1zu; j < n_columns; ++j) {
            const auto dlogx = std::log(SB_k_over_gm1[j] / SB_k_over_gm1[j - 1]);
            const auto dsigma = 0.5 * (kdsdk[j - 1] + kdsdk[j]) * dlogx;
            sigma_ttl += dsigma;
            data[j + 1][i] = sigma_ttl;
         }

         // Normalize the CDF.  This makes data[1] = 0 and data[n_columns] = 1.
         // The total cross section used by the collision probability is made
         // self-consistent with this same log-integrated table.
         if (sigma_ttl > 0.0) {
            for (auto j = 0zu; j < n_columns; ++j) {
               data[j + 1][i] /= sigma_ttl;
            }
         }

         data[n_columns + 1][i] = sigma_ttl * 1.0e-28; // barns -> m^2
      }
      // for (auto i = 0zu; i < n_energies; ++i) {
      //    for (auto j = 0zu; j < n_columns + 2; ++j) {
      //       std::print("{}, ", data[j][i]);
      //    }
      //    std::println();
      // }
   }

   bool is_outofbounds(const auto e) const {
      return e <= data[0][0] or e > data[0][data[0].size() - 1];
   }

   auto lerp(const auto energy, const auto U) const {
      const auto n_columns = SB_k_over_gm1.size();

      auto e_begin = data[0].begin();
      auto e_end = data[0].end();
      auto e_lb_index = std::lower_bound(e_begin, e_end, energy) - e_begin;

      auto de_lb = std::abs(data[0][static_cast<size_t>(e_lb_index)] - energy);
      auto de_lbm1 = std::abs(data[0][static_cast<size_t>(e_lb_index - 1)] - energy);

      const auto e_nearest_index = static_cast<size_t>((de_lb < de_lbm1) ? e_lb_index : (e_lb_index - 1));

      // Clamp the random variate into the open unit interval so that roundoff
      // never asks for a value below the first CDF entry or above the last.
      const auto Uc = std::min(std::max(static_cast<double>(U), 0.0), 1.0 - 1.0e-14);

      // data[1] is CDF=0 at x_min, and data[n_columns] is CDF=1 at x_max.
      // Find the first CDF point that exceeds the random variate.
      size_t k_column = 1;
      auto sigma_cdf = data[k_column + 1][e_nearest_index];
      while (sigma_cdf < Uc and k_column < n_columns - 1) {
         ++k_column;
         sigma_cdf = data[k_column + 1][e_nearest_index];
      }

      const auto sigma_cdf_m1 = data[k_column][e_nearest_index];
      const auto denom = sigma_cdf - sigma_cdf_m1;
      const auto frac = denom > 0.0 ? (Uc - sigma_cdf_m1) / denom : 0.0;

      // Interpolate in log photon-energy fraction, not linearly in x.  This is
      // appropriate because the table is logarithmically spaced and represents
      // k*d(sigma)/dk = d(sigma)/dln(k).
      const auto logx0 = std::log(SB_k_over_gm1[k_column - 1]);
      const auto logx1 = std::log(SB_k_over_gm1[k_column]);
      const auto logx = logx0 + frac * (logx1 - logx0);

      return std::exp(logx);
   }

   auto lerp_cumulative(const auto e) const {
      const auto n_col = data.size() - 1;
      const auto n_row = data[0].size() - 1;
      if (e < data[0][1]) { return data[n_col][0]; }
      if (e > data[0][n_row]) { return data[n_col][n_row]; }
      const auto upper = std::ranges::upper_bound(data[0], e);
      const auto idx = std::ranges::distance(data[0].cbegin(), upper) - 1;
      const auto slope = (e - data[0][idx]) / (data[0][idx + 1] - data[0][idx]);
      return std::lerp(data[n_col][idx], data[n_col][idx + 1], slope);
   }
}; // end struct MultiTable



} // end namespace tf::interp
#endif //INTERPOLATION_HPP
