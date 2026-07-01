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
#include <boost/math/special_functions/lambert_w.hpp>

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

struct BremTable {
   std::vector<std::vector<double>> data;
   std::array<double, 14> SB_k_over_gm1 = {1.0e-7, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.97, 0.99, 1.0};

   explicit BremTable() = default;
   explicit BremTable(const std::string& filepath)
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
         // tabulated values of x = k/(gamma-1).  The photon-number cross-section
         // is therefore the integral of k*d(sigma)/dk over dln(k),
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
         auto sigma_ttl = 0.0;

         for (auto j = 1zu; j < n_columns; ++j) {
            const auto dlogx = std::log(SB_k_over_gm1[j] / SB_k_over_gm1[j - 1]);
            const auto dsigma = 0.5 * (kdsdk[j - 1] + kdsdk[j]) * dlogx;
            sigma_ttl += dsigma;
            data[j + 1][i] = sigma_ttl;
         }

         // Normalize the CDF.  This makes data[1] = 0 and data[n_columns] = 1.
         // The total cross-section used by the collision probability is made
         // self-consistent with this same log-integrated table.
         if (sigma_ttl > 0.0) {
            for (auto j = 0zu; j < n_columns; ++j) {
               data[j + 1][i] /= sigma_ttl;
            }
         }

         data[n_columns + 1][i] = sigma_ttl * 1.0e-28; // barns -> m^2
      }
   }

   auto is_outofbounds(const auto e) const -> bool {
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

struct BremTFDTable {
   static constexpr auto num_electron_energies = 10zu; // rows
   static constexpr auto num_photon_energies = 10zu; // columns
   static constexpr auto k_over_gm1_cutoff = 1.0e-7;
   static constexpr auto delta_debye_cutoff = 0.1;
   static constexpr std::array<double, 10> gl_roots = {-0.97390653, -0.86506337, -0.67940957, -0.43339539, -0.14887434,
                                                        0.14887434, 0.43339539, 0.67940957, 0.86506337, 0.97390653};
   static constexpr std::array<double, 10> gl_weights = {0.06667134, 0.14945135, 0.21908636, 0.26926672, 0.29552422,
                                                         0.29552422, 0.26926672, 0.21908636, 0.14945135, 0.06667134};

   struct CellCDF {
      CellCDF()
      : cdf(num_electron_energies, num_photon_energies),
        total_cross_sections(num_electron_energies),
        prev_debye(0.0)
      {}

      Array2D<double> cdf;
      std::vector<double> total_cross_sections{};
      double prev_debye;
   };

   std::vector<double> energies{};
   std::vector<double> k_over_gm1{};

   std::unordered_map<std::size_t, CellCDF> cdfs{};

   double atomic_number{};
   double eta_D{};
   double eta_TF{};
   double lTF{};

   explicit BremTFDTable() = default;

   BremTFDTable(const auto Z, const auto charge, const auto Emin, const auto Emax)
   : energies(math::logspace(std::log10(Emin), std::log10(Emax), num_electron_energies)),
     k_over_gm1(math::logspace(std::log10(k_over_gm1_cutoff), std::log10(1.0 - k_over_gm1_cutoff), num_photon_energies)),
     atomic_number(static_cast<double>(Z)),
     eta_D(charge / atomic_number),
     eta_TF(1.0 - eta_D),
     lTF(4.0 * constants::pi * constants::eps0 * math::SQR(constants::h_bar) / (constants::m_e * math::SQR(constants::q_e) * std::pow(atomic_number, 1.0 / 3.0)))
   {}

   auto ElwertFactor(const auto k, const auto gamma) const -> double {
      const auto b1 = std::sqrt(1.0 - 1.0 / math::SQR(gamma));
      const auto b2 = std::sqrt(1.0 - 1.0 / math::SQR(gamma - k));

      const auto cond = atomic_number * constants::alpha_fine * (1.0 / b2 - 1.0 / b1);
      if (cond < 100.0) {
         return (b1 / b2) * (1.0 - std::exp(-2.0 * constants::pi * atomic_number * constants::alpha_fine / b1))
                          / (1.0 - std::exp(-2.0 * constants::pi * atomic_number * constants::alpha_fine / b2));
      }
      return 1.0;
   }

   auto compute_nonrelativistic_dsdk(const auto k, const auto gamma, const auto lD) const -> double {
      const auto p1 = std::sqrt(math::SQR(gamma) - 1.0);
      const auto p2 = std::sqrt(math::SQR(gamma - k) - 1.0);

      const auto dp = p1 + p2;
      const auto dm = p1 - p2;

      const auto lD2 = math::SQR(lD);
      const auto lTF2 = math::SQR(lTF);

      const auto dpl2_TF = 1.0 / (math::SQR(dp) * lTF2 + 1.0);
      const auto dml2_TF = 1.0 / (math::SQR(dm) * lTF2 + 1.0);
      const auto T_TF = 0.5 * math::SQR(eta_TF) * (std::log(dml2_TF / dpl2_TF) + dpl2_TF - dml2_TF);

      const auto dpl2_D = 1.0 / (math::SQR(dp) * lD2 + 1.0);
      const auto dml2_D = 1.0 / (math::SQR(dm) * lD2 + 1.0);
      const auto T_D = 0.5 * math::SQR(eta_D) * (std::log(dml2_D / dpl2_D) + dpl2_D - dml2_D);

      const auto elwert = ElwertFactor(k, gamma - 1.0);
      const auto T_coef = 16.0 * math::SQR(atomic_number * constants::r_e) * constants::alpha_fine / (3.0 * k * math::SQR(p1));
      const auto T_coupling = eta_TF * eta_D / (lD2 - lTF2) * (lTF2 * std::log(dpl2_D / dml2_D) + lD2 * std::log(dml2_TF / dpl2_TF));

      return elwert * T_coef * (T_TF + T_D + T_coupling);
   }

   auto I1_screening(const auto d, const auto l) const -> double {
      const auto l2 = math::SQR(l);
      const auto ld = l * d;
      const auto T1 = ld * (std::atan(ld) - std::atan(l));
      const auto T2 = -0.5 * math::SQR(l * (1.0 - d)) / (1.0 + l2);
      const auto T3 = 0.5 * std::log((1.0 + l2) / (1.0 + math::SQR(ld)));
      return T1 + T2 + T3;
   }

   auto I2_screening(const auto d, const auto l) const -> double {
      const auto l2 = math::SQR(l);
      const auto ld = l * d;
      const auto ld2 = math::SQR(ld);
      const auto T1 = 4.0 * math::CUBE(ld) * (std::atan(ld) - std::atan(l));
      const auto T2 = (1.0 + 3.0 * ld2) * std::log((1.0 + l2) / (1.0 + ld2));
      const auto T3 = (6.0 * l2 * ld2) * std::log(d) / (1.0 + l2);
      const auto T4 = l2 * (d - 1.0) * (d + 1.0 - 4.0 * ld2) / (1.0 + l2);
      return 0.5 * (T1 + T2 + T3 + T4);
   }

   auto compute_relativistic_dsdk(const auto k, const auto gamma, const auto lR) const -> double {
      const auto gmk = gamma - k;
      const auto d = k / (2.0 * gamma * gmk);
      const auto g_ratio = gmk / gamma;
      const auto T_coef = 4.0 * math::SQR(atomic_number * constants::r_e) * constants::alpha_fine / k;
      const auto T1 = (1.0 + math::SQR(g_ratio)) * (I1_screening(d, lR) + 1.0);
      const auto T2 = (2.0 / 3.0) * g_ratio * (I2_screening(d, lR) + 5.0 / 6.0);
      return T_coef * (T1 - T2);
   }

   [[nodiscard]] auto zetaRiemann(const auto s, const auto n) const -> double {
      auto zeta = 0.0;
      for (auto i = 0zu; i < n; ++i) {
         zeta += 1.0 / std::pow(static_cast<double>(i) + 1.0, s);
      }
      return zeta;
   }

   auto compute_ultrarelativistic_dsdk(const auto k, const auto gamma, const auto lR) const -> double {
      const auto gmk = gamma - k;
      const auto d = k / (2.0 * gamma * gmk);
      const auto g_ratio = gmk / gamma;

      const auto aZ2 =  math::SQR(constants::alpha_fine * atomic_number);
      const auto t1 = aZ2 / (1.0 + aZ2);

      auto t2 = 0.0;
      for (auto i = 1zu; i < 5; ++i) {
         // const auto zeta = zetaRiemann(2 * i + 1, 10);
         const auto zeta = std::riemann_zeta(static_cast<double>(2 * i + 1));
         // assert(std::abs(zeta - zeta2) < 1.0e-8);
         t2 += std::pow(-aZ2, i) * (zeta - 1.0);
      }
      const auto f_coulomb = t1 * t2;

      const auto T_coef = 4.0 * math::SQR(atomic_number * constants::r_e) * constants::alpha_fine / k;
      const auto T1 = (1.0 + math::SQR(g_ratio)) * (I1_screening(d, lR) + 1.0 - f_coulomb);
      const auto T2 = (2.0 / 3.0) * g_ratio * (I2_screening(d, lR) + 5.0 / 6.0 - f_coulomb);

      return T_coef * (T1 - T2);
   }

   auto compute_dsdk(const auto k, const auto gm1, const auto lD, const auto lR) const -> double {
      const auto gamma = gm1 + 1.0;
      assert(gamma >= 1.0);

      if (k < 0.0 or k > gm1) { return 0.0; }

      if (gamma <= 2.0) {
         return compute_nonrelativistic_dsdk(k, gamma, lD);
      }
      if (gamma > 2.0 and gamma <= 100.0) {
         return compute_relativistic_dsdk(k, gamma, lR);
      }
      if (gamma > 100.0) {
         return compute_ultrarelativistic_dsdk(k, gamma, lR);
      }
      return 0.0;
   }

   auto integrate_total_cs(const auto gm1, const auto lD, const auto lR) const -> double {
      constexpr auto degree = 10zu;

      const auto a = std::log(k_over_gm1_cutoff * gm1);
      const auto b = std::log(gm1);

      auto ttl_cs = 0.0;
      for (auto i = 0zu; i < degree; ++i) {
         const auto k = std::exp(0.5 * (gl_roots[i] + 1.0) * (b - a) + a);
         const auto k_dsdk = k * compute_dsdk(k, gm1, lD, lR);
         ttl_cs += gl_weights[i] * k_dsdk * 0.5 * (b - a);
      }
      return ttl_cs;
   }

   auto reducedPotentialLength(const auto lD) const -> double {
      const auto lTF2 = math::SQR(lTF);
      const auto lD2 = math::SQR(lD);
      const auto T1 = 0.5 * math::SQR(eta_TF) * ((1.0 + lTF2) * std::log(1.0 + lTF2) - lTF2) / (1.0 + lTF2);
      auto T2 = 0.0;
      auto T3 = 0.0;

      if (lD > 0.0) { // todo: does this every NOT trigger? can lD every be 0?
         T2 = 0.5 * math::SQR(eta_D) * ((1.0 + lD2) * std::log(1.0 + lD2) - lD2) / (1.0 + lD2);
         T3 = eta_TF * eta_D * (lD2 * std::log(1.0 + lTF2) - lTF2 * std::log(1.0 + lD2)) / (lD2 - lTF2);
      }

      const auto I_TFD = T1 + T2 + T3;
      const auto x = 1.0 + 2.0 * I_TFD;
      const auto w = boost::math::lambert_w0(-1.0 * std::exp(-x));
      return std::sqrt(std::exp(w + x) - 1.0);
   }

   auto computeCDF(auto& cell, const auto l_debye) {
      const auto lD = l_debye / constants::r_c;
      const auto lR = reducedPotentialLength(lD);

      for (auto i = 0zu; i < num_electron_energies; ++i) {
         const auto gm1 = energies[i] * constants::q_e / constants::m_e_c_sqr;

         auto k_jm1 = 0.0;
         auto sigma_ttl = 0.0;

         for (auto j = 0zu; j < num_photon_energies; ++j) {
            const auto k = k_over_gm1[j] * gm1;
            const auto dk = k - k_jm1;
            k_jm1 = k;

            const auto sigma = dk * compute_dsdk(k, gm1, lD, lR);
            cell.cdf(i, j) = sigma;
            sigma_ttl += sigma;
         }

         cell.total_cross_sections[i] = integrate_total_cs(gm1, lD, lR);
         cell.cdf(i, 0) /= sigma_ttl;
         for (auto j = 1zu; j < num_photon_energies; ++j) {
            cell.cdf(i, j) = cell.cdf(i, j - 1) + cell.cdf(i, j) / sigma_ttl;
         }
      }
   }

   void updateCDF(const auto cid, const auto l_debye) {
      if (not cdfs.contains(cid)) {
         cdfs.emplace(cid, CellCDF{});
      }

      auto& cdf = cdfs.at(cid);
      const auto delta_debye = std::abs(cdf.prev_debye - l_debye) / cdf.prev_debye;

      if (delta_debye > delta_debye_cutoff) {
         cdf.prev_debye = l_debye;
         computeCDF(cdf, l_debye);
      }
   }

   auto is_outofbounds(const auto e) const -> bool {
      return e <= energies[0] or e > energies[energies.size() - 1];
   }

   auto interpolate(const auto e, const auto cid) const {
      const auto& cell = cdfs.at(cid);
      const auto n_row = energies.size() - 1;
      if (e < energies[0]) { return cell.total_cross_sections[0]; }
      if (e > energies[n_row]) { return cell.total_cross_sections[n_row]; }
      const auto upper = std::ranges::upper_bound(energies, e);
      const auto idx = std::ranges::distance(energies.cbegin(), upper) - 1;
      const auto slope = (e - energies[idx]) / (energies[idx + 1] - energies[idx]);
      return std::lerp(cell.total_cross_sections[idx], cell.total_cross_sections[idx + 1], slope);
   }

   auto lerp(const auto e, const auto U, const auto cid) const {
      const auto& cell = cdfs.at(cid);

      const auto n_columns = k_over_gm1.size();

      const auto upper = std::ranges::upper_bound(energies, e);
      const auto idx = std::ranges::distance(energies.cbegin(), upper) - 1;

      size_t k_column = 1;
      auto sigma_cdf = cell.cdf(k_column + 1, idx);
      while (sigma_cdf < U and k_column < n_columns - 1) {
         ++k_column;
         sigma_cdf =  cell.cdf(k_column + 1, idx);
      }

      const auto sigma_cdf_m1 =  cell.cdf(k_column, idx);
      const auto denom = sigma_cdf - sigma_cdf_m1;
      const auto frac = denom > 0.0 ? (U - sigma_cdf_m1) / denom : 0.0;

      const auto logx0 = std::log(k_over_gm1[k_column - 1]);
      const auto logx1 = std::log(k_over_gm1[k_column]);
      const auto logx = logx0 + frac * (logx1 - logx0);

      return std::exp(logx);
   }
};

} // end namespace tf::interp
#endif //INTERPOLATION_HPP
