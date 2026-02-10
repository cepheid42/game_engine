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

namespace tf::interp
{
/*---------------------------------------------------------------/
/-                     Interpolation Classes                    -/
/---------------------------------------------------------------*/
template<int D>
constexpr auto rotateOrigin(const auto x, const auto y, const auto z) {
   if      constexpr (D == 0) { return std::array{y, z, x}; } // 1, 2, 0
   else if constexpr (D == 1) { return std::array{z, x, y}; } // 2, 0, 1
   else                       { return std::array{x, y, z}; } // 0, 1, 2
}

struct NGP {
   static constexpr int         Begin   = -1;
   static constexpr int         End     = -1;
   static constexpr std::size_t Order   = 0;
   static constexpr std::size_t Support = 1;

   static constexpr auto eval(const double) {
      return 1.0;
   }

   static constexpr auto shape_array(const double) {
      return std::array{1.0};
   }

   static constexpr auto ds_array(const auto, const auto&) {
      return std::array{1.0};
   }
};

struct CIC {
   static constexpr int         Begin   = -1;
   static constexpr int         End     = 0;
   static constexpr std::size_t Order   = 1;
   static constexpr std::size_t Support = 2;

   static constexpr auto eval(const double x) {
      return 1.0 - std::abs(x);
   }

   static constexpr auto shape_array(const double x) {
      return std::array{eval(x - Begin), eval(x)};
   }

   static constexpr auto ds_array(const auto x1, const auto& s0) {
      return std::array{eval(x1 - Begin) - s0[0],
                        eval(x1)         - s0[1]};
   }
};

struct TSC {
   static constexpr int         Begin   = -1;
   static constexpr int         End     = 1;
   static constexpr std::size_t Order   = 2;
   static constexpr std::size_t Support = 3;

   static constexpr auto innerRadius(const auto x) {
      return 0.75 - math::SQR(x);
   }

   static constexpr auto outerRadius(const auto x) {
      return 0.5 * math::SQR(1.5 - x);
   }

   static constexpr auto eval(const auto x) {
      const auto absx = std::abs(x);
      return absx <= 0.5 ? innerRadius(absx) : outerRadius(absx);
   }

   static constexpr auto shape_array(const auto x) {
      return std::array{eval(x - Begin), eval(x), eval(x - End)};
   }

   static constexpr auto ds_array(const auto x1, const auto& s0) {
      return std::array{eval(x1 - Begin) - s0[0],
                        eval(x1)         - s0[1],
                        eval(x1 - End)   - s0[2]};
   }
};

struct PQS {
   static constexpr int         Begin   = -1;
   static constexpr int         End     = 2;
   static constexpr std::size_t Order   = 3;
   static constexpr std::size_t Support = 4;

   static constexpr auto innerRadius(const double x) {
      return (2.0 / 3.0) - math::SQR(x) + 0.5 * math::CUBE(x);
   }

   static constexpr auto outerRadius(const double x) {
      return (1.0 / 6.0) * math::CUBE(2.0 - x);
   }

   static constexpr auto eval(const double x) {
      const auto absx = std::abs(x);
      return absx <= 1.0 ? innerRadius(x) : outerRadius(absx);
   }

   static constexpr auto shape_array(const double x) {
      return std::array{eval(x - Begin), eval(x), eval(x - 1.0), eval(x - End)};
   }

   static constexpr auto ds_array(const auto x1, const auto& s0) {
      return std::array{eval(x1 - Begin) - s0[0],
                        eval(x1)         - s0[1],
                        eval(x1 - 1.0)   - s0[2],
                        eval(x1 - End)   - s0[2]};
   }
};

template<int ShapeOrder> struct InterpolationShape;
template<> struct InterpolationShape<0> { using Type = NGP; };
template<> struct InterpolationShape<1> { using Type = CIC; };
template<> struct InterpolationShape<2> { using Type = TSC; };
template<> struct InterpolationShape<3> { using Type = PQS; };

template<typename Outer, typename Middle, typename Inner>
struct InterpolationStrategy {
   using OuterShape = Outer;
   using MiddleShape = Middle;
   using InnerShape = Inner;
};


/*---------------------------------------------------------------/
/-                     Table Interpolation                      -/
/---------------------------------------------------------------*/
struct Table {
   std::vector<double> xs{};
   std::vector<double> ys{};

   explicit Table() = default;

   explicit Table(const std::string& filepath, const auto xscale=1.0, const auto yscale=1.0) {
      std::ifstream file(std::string{sim_path} + filepath);

      if (!file.is_open()) {
         throw std::runtime_error("Unable to open cross section table: " + filepath);
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
      std::ifstream file(std::string{sim_path} + filepath);

      if (!file.is_open()) {
         throw std::runtime_error("Unable to open cross section table: " + filepath);
      }

      double val{};
      std::string line;
      while (std::getline(file, line)) {
         std::istringstream buffer(line);
         std::size_t col{0zu};
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

         const auto energy_eV = data[0][i];
         const auto v = constants::c * std::sqrt(1.0 - 1.0 / math::SQR(1.0 + constants::q_e * energy_eV / constants::m_e_c_sqr));
         const auto gm1 = 1.0 / std::sqrt(1.0 - math::SQR(v / constants::c)) - 1.0;

         const auto sigma_ttl = data[n_columns + 1][i];

         for (auto j = 1zu; j < n_columns; ++j) {
            const auto k = SB_k_over_gm1[j - 1] * gm1;
            double dk;
            if (j == 1) {
               dk = k;
            } else {
               dk = SB_k_over_gm1[j] * gm1 - k;
            }

            auto sigma_norm = data[j][i] * dk / k / sigma_ttl;

            if (j == 1) {
               data[j][i] = sigma_norm;
            } else {
               data[j][i] = data[j - 1][i] + sigma_norm;
            }
         }
         data[n_columns + 1][i] *= 1.0e-28; // barns -> m^2
      }
      // for (auto i = 0zu; i < n_energies; ++i) {
      //    for (auto j = 0zu; j < n_columns + 2; ++j) {
      //       std::print("{}, ", data[j][i]);
      //    }
      //    std::println();
      // }
   }

   auto lerp(const auto energy, const auto U) const {
      const auto n_columns = SB_k_over_gm1.size();

      auto e_begin = data[0].begin();
      auto e_end = data[0].end();
      auto e_lb_index = std::lower_bound(e_begin, e_end, energy) - e_begin;

      auto de_lb = std::abs(data[0][static_cast<size_t>(e_lb_index)] - energy);
      auto de_lbm1 = std::abs(data[0][static_cast<size_t>(e_lb_index - 1)] - energy);

      auto e_nearest_index = (de_lb < de_lbm1) ? e_lb_index : (e_lb_index - 1);

      size_t k_column = 0; // Starts at one since energy is in first column of table
      auto sigma_cdf = data[k_column + 1][static_cast<size_t>(e_nearest_index)];
      while (sigma_cdf < U and k_column < n_columns - 1) {
         ++k_column;
         sigma_cdf = data[k_column + 1][static_cast<size_t>(e_nearest_index)];
      } // endwhile()

      const auto k_over_gm1_i = SB_k_over_gm1[k_column];

      double k_over_gm1{};
      if (k_column == 0) {
         k_over_gm1 = k_over_gm1_i * U / sigma_cdf;
      } else {
         const auto k_over_gm1_m1 = SB_k_over_gm1[k_column - 1];
         const auto sigma_norm_m1 = data[k_column][static_cast<size_t>(e_nearest_index)];
         const auto frac = (U - sigma_norm_m1) / (sigma_cdf - sigma_norm_m1);
         k_over_gm1 = k_over_gm1_m1 + (k_over_gm1_i - k_over_gm1_m1) * frac;
      }

      return k_over_gm1;
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
