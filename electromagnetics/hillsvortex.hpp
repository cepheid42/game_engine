#ifndef HILLSVORTEX_HPP
#define HILLSVORTEX_HPP

#include "program_params.hpp"
#include "math_utils.hpp"

#include <cmath>

namespace tf::electromagnetics {

struct HillsFieldGenerator {
   static void fillEMData(auto& emdata, const auto B0, const auto rs, const auto kk) {
      using math::SQR;

      for (size_t i = 0; i < emdata.Bx_app.nx(); i++) {
         for (size_t j = 0; j < emdata.Bx_app.ny(); j++) {
            for (size_t k = 0; k < emdata.Bx_app.nz(); k++) {
               const auto xx = x_range[0] + (static_cast<double>(i) + 0.5) * dx;
               const auto yy = y_range[0] + (static_cast<double>(j) + 0.0) * dy;
               const auto zz = z_range[0] + (static_cast<double>(k) + 0.0) * dz;
               const auto rr = std::sqrt(SQR(xx) + SQR(zz));
               const auto theta = std::atan2(zz, xx);
               emdata.Bx_app(i, j, k) = getHillsBr(rr, yy, B0, rs, kk) * std::cos(theta);
            }
         }
      }

      for (size_t i = 0; i < emdata.By_app.nx(); i++) {
         for (size_t j = 0; j < emdata.By_app.ny(); j++) {
            for (size_t k = 0; k < emdata.By_app.nz(); k++) {
               const auto xx = x_range[0] + (static_cast<double>(i) + 0.0) * dx;
               const auto yy = y_range[0] + (static_cast<double>(j) + 0.5) * dy;
               const auto zz = z_range[0] + (static_cast<double>(k) + 0.0) * dz;
               const auto rr = std::sqrt(SQR(xx) + SQR(zz));
               emdata.By_app(i, j, k) = getHillsBz(rr, yy, B0, rs, kk);
            }
         }
      }

      for (size_t i = 0; i < emdata.Bz_app.nx(); i++) {
         for (size_t j = 0; j < emdata.Bz_app.ny(); j++) {
            for (size_t k = 0; k < emdata.Bz_app.nz(); k++) {
               const auto xx = x_range[0] + (static_cast<double>(i) + 0.0) * dx;
               const auto yy = y_range[0] + (static_cast<double>(j) + 0.0) * dy;
               const auto zz = z_range[0] + (static_cast<double>(k) + 0.5)* dz;
               const auto rr = std::sqrt(SQR(xx) + SQR(zz));
               const auto theta = std::atan2(zz, xx);
               emdata.Bz_app(i, j, k) = getHillsBr(rr, yy, B0, rs, kk) * std::sin(theta);
            }
         }
      }
   }


   static auto getHillsBz(const auto r, auto z, const auto B0, const auto a, const auto kk) {
      using math::SQR;
      using math::CUBE;
      z = z / kk;
      const auto rsqr = SQR(r);
      const auto zsqr = SQR(z);
      const auto asqr = SQR(a);
      const auto acbe = CUBE(a);
      const auto R = std::sqrt(rsqr + zsqr);
      return R >= a ? B0 * (1.0 - acbe / std::pow(R, 1.5) + 1.5 * acbe * rsqr / std::pow(R, 2.5))
                    : B0 * 1.5 / asqr * (-1.0 * asqr + 2.0 * rsqr + zsqr);
   }

   static auto getHillsBr(const auto r, auto z, const auto B0, const auto a, const auto kk) {
      using math::SQR;
      using math::CUBE;
      z = z / kk;
      const auto rsqr = SQR(r);
      const auto zsqr = SQR(z);
      const auto asqr = SQR(a);
      const auto acbe = CUBE(a);
      const auto R = std::sqrt(rsqr + zsqr);
      return R >= a ? -1.5 * B0 * acbe * r * z / std::pow(R, 2.5)
                    : -1.0 * B0 * 1.5 / asqr * r * z;
   }
};

}


#endif //HILLSVORTEX_HPP
