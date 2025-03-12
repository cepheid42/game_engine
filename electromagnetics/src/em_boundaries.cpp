// #include <print>
// #include <vector>
//
// #include "em_boundaries.hpp"
//
// #include "program_params.hpp"
// #include "constants.hpp"
// #include "math_utils.hpp"
//
// namespace tf::electromagnetics {
//   void PMLData::init_coefficients() {
//     constexpr auto coef1 = -dt / constants::eps0;
//
//     auto hstep = 1.0 / (2.0 * static_cast<double>(PMLDepth));
//     std::vector<double> d = math::linspace(1.0, 0.0, PMLDepth, false);
//     std::vector<double> dh = math::linspace(1.0 - hstep, -hstep, PMLDepth, false);
//
//     for (auto e: d) {
//       std::print("{}, ", e);
//     }
//     std::println();
//     for (auto e: dh) {
//       std::print("{}, ", e);
//     }
//     std::println();
//
//     exit(0);
//   }
//
// }
