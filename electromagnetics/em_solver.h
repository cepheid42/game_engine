//
// Created by cepheid on 6/28/24.
//

#ifndef EM_SOLVER_H
#define EM_SOLVER_H

#include <tfsf.h>

#include "aydenstuff/array.h"
// #include "em_updates.h"
// #include "em_sources.h"


namespace tf::electromagnetics
{
  template<typename EXI, typename EYI, typename EZI,
           typename HXI, typename HYI, typename HZI,
           typename BCX0, typename BCX1,
           typename BCY0, typename BCY1,
           typename BCZ0, typename BCZ1>
  struct Electromagnetics {
    using value_t = typename EXI::value_t;
    using dimension_t = typename EXI::dimension_t;
    using empty_t = tf::types::EmptyArray<value_t, dimension_t::value>;

    static constexpr empty_t empty{};

    static void updateE(auto& emdata) {
      EXI::apply(emdata.Ex, emdata.Hz, emdata.Hy, emdata.Jx, emdata.Cexe, emdata.Cexhz, emdata.Cexhy, emdata.Cjx, {0, 0, 1, 1, 1, 1});
      EYI::apply(emdata.Ey, emdata.Hx, emdata.Hz, emdata.Jy, emdata.Ceye, emdata.Ceyhx, emdata.Ceyhz, emdata.Cjy, {1, 1, 0, 0, 1, 1});
      EZI::apply(emdata.Ez, emdata.Hy, emdata.Hx, emdata.Jz, emdata.Ceze, emdata.Cezhy, emdata.Cezhx, emdata.Cjz, {1, 1, 1, 1, 0, 0});
    }

    static void updateH(auto& emdata) {
      HXI::apply(emdata.Hx, emdata.Ey, emdata.Ez, empty, emdata.Chxh, emdata.Chxey, emdata.Chxez, empty, {0, 0, 0, 0, 0, 0});
      HYI::apply(emdata.Hy, emdata.Ez, emdata.Ex, empty, emdata.Chyh, emdata.Chyez, emdata.Chyex, empty, {0, 0, 0, 0, 0, 0});
      HZI::apply(emdata.Hz, emdata.Ex, emdata.Ey, empty, emdata.Chzh, emdata.Chzex, emdata.Chzey, empty, {0, 0, 0, 0, 0, 0});
    }

    static void updateE_bcs(auto& emdata, auto& bcdata) {
      BCX0::Ex::updateE(bcdata.x0.Ex, emdata.Ex, empty, empty);
      BCX0::Ey::updateE(bcdata.x0.Ey, emdata.Ey, emdata.Hz, emdata.Ceyhz);
      BCX0::Ez::updateE(bcdata.x0.Ez, emdata.Ez, emdata.Hy, emdata.Cezhy);

      BCX1::Ex::updateE(bcdata.x1.Ex, emdata.Ex, empty, empty);
      BCX1::Ey::updateE(bcdata.x1.Ey, emdata.Ey, emdata.Hz, emdata.Ceyhz);
      BCX1::Ez::updateE(bcdata.x1.Ez, emdata.Ez, emdata.Hy, emdata.Cezhy);

      BCY0::Ex::updateE(bcdata.y0.Ex, emdata.Ex, emdata.Hz, emdata.Cexhz);
      BCY0::Ey::updateE(bcdata.y0.Ey, emdata.Ey, empty, empty);
      BCY0::Ez::updateE(bcdata.y0.Ez, emdata.Ez, emdata.Hx, emdata.Cezhx);

      BCY1::Ex::updateE(bcdata.y1.Ex, emdata.Ex, emdata.Hz, emdata.Cexhz);
      BCY1::Ey::updateE(bcdata.y1.Ey, emdata.Ey, empty, empty);
      BCY1::Ez::updateE(bcdata.y1.Ez, emdata.Ez, emdata.Hx, emdata.Cezhx);

      BCZ0::Ex::updateE(bcdata.z0.Ex, emdata.Ex, emdata.Hy, emdata.Cexhy);
      BCZ0::Ey::updateE(bcdata.z0.Ey, emdata.Ey, emdata.Hx, emdata.Ceyhx);
      BCZ0::Ez::updateE(bcdata.z0.Ez, emdata.Ez, empty, empty);

      BCZ1::Ex::updateE(bcdata.z1.Ex, emdata.Ex, emdata.Hy, emdata.Cexhy);
      BCZ1::Ey::updateE(bcdata.z1.Ey, emdata.Ey, emdata.Hx, emdata.Ceyhx);
      BCZ1::Ez::updateE(bcdata.z1.Ez, emdata.Ez, empty, empty);
    }

    static void updateH_bcs(auto& emdata, auto& bcdata) {
      BCX0::Hx::updateH(bcdata.x0.Hx, emdata.Hx, empty, empty);
      BCX0::Hy::updateH(bcdata.x0.Hy, emdata.Hy, emdata.Ez, emdata.Chyez);
      BCX0::Hz::updateH(bcdata.x0.Hz, emdata.Hz, emdata.Ey, emdata.Chzey);

      BCX1::Hx::updateH(bcdata.x1.Hx, emdata.Hx, empty, empty);
      BCX1::Hy::updateH(bcdata.x1.Hy, emdata.Hy, emdata.Ez, emdata.Chyez);
      BCX1::Hz::updateH(bcdata.x1.Hz, emdata.Hz, emdata.Ey, emdata.Chzey);

      BCY0::Hx::updateH(bcdata.y0.Hx, emdata.Hx, emdata.Ez, emdata.Chxez);
      BCY0::Hy::updateH(bcdata.y0.Hy, emdata.Hy, empty, empty);
      BCY0::Hz::updateH(bcdata.y0.Hz, emdata.Hz, emdata.Ex, emdata.Chzex);

      BCY1::Hx::updateH(bcdata.y1.Hx, emdata.Hx, emdata.Ez, emdata.Chxez);
      BCY1::Hy::updateH(bcdata.y1.Hy, emdata.Hy, empty, empty);
      BCY1::Hz::updateH(bcdata.y1.Hz, emdata.Hz, emdata.Ex, emdata.Chzex);

      BCZ0::Hx::updateH(bcdata.z0.Hx, emdata.Hx, emdata.Ey, emdata.Chxey);
      BCZ0::Hy::updateH(bcdata.z0.Hy, emdata.Hy, emdata.Ex, emdata.Chyex);
      BCZ0::Hz::updateH(bcdata.z0.Hz, emdata.Hz, empty, empty);

      BCZ1::Hx::updateH(bcdata.z1.Hx, emdata.Hx, emdata.Ey, emdata.Chxey);
      BCZ1::Hy::updateH(bcdata.z1.Hy, emdata.Hy, emdata.Ex, emdata.Chyex);
      BCZ1::Hz::updateH(bcdata.z1.Hz, emdata.Hz, empty, empty);
    }


    static void advance(auto q, auto& emdata, auto& bcdata, auto& src) {
      updateH(emdata);
      updateH_bcs(emdata, bcdata);

      src.apply(emdata, q);
      // emdata.Ez(60, 60) += tf::electromagnetics::sources::detail::ricker(q);

      updateE(emdata);
      updateE_bcs(emdata, bcdata);

      updateH(emdata);
      updateH_bcs(emdata, bcdata);
    }
  };
} // end namespace tf::electromagnetics
#endif //EM_SOLVER_H
