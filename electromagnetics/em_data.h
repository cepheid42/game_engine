//
// Created by cepheid on 9/23/24.
//

#ifndef EM_DATA_H
#define EM_DATA_H

#include "core/typelist.h"
#include "aydenstuff/array.h"
#include "core/debug.h"
#include "em_emtpyarray.h"
#include "em_traits.h"

using tf::types::Array1D;
using tf::types::Array2D;

template<FieldComponent EXF, FieldComponent EYF, FieldComponent EZF, FieldComponent HXF, FieldComponent HYF, FieldComponent HZF>
struct EMData {
  using value_t = typename EXF::arr_t::value_t;
  using dimension_t = typename EXF::arr_t::dimension_t;
  
  using empty_t = EmptyArray<value_t, dimension_t::value>;
  
  using ex_t = typename EXF::arr_t;
  using ey_t = typename EYF::arr_t;
  using ez_t = typename EZF::arr_t;
  using hx_t = typename HXF::arr_t;
  using hy_t = typename HYF::arr_t;
  using hz_t = typename HZF::arr_t;
  
  EMData() = default;
  
  explicit EMData(size_t nx)
  : Ex{nx}, Jx{nx}, Cexe{nx}, Cexh{nx}, Cjx{nx},
    Ey{nx}, Jy{nx}, Ceye{nx}, Ceyh{nx}, Cjy{nx},
    Ez{nx}, Jz{nx}, Ceze{nx}, Cezh{nx}, Cjz{nx},
    Hx{nx}, Chxe{nx}, Chxh{nx},
    Hy{nx}, Chye{nx}, Chyh{nx},
    Hz{nx}, Chze{nx}, Chzh{nx}
  {
    DBG("EMData::EMData(nx)", nx);
  }
  
  explicit EMData(size_t nx, size_t ny)
  : Ex{nx, ny}, Jx{nx, ny}, Cexe{nx, ny}, Cexh{nx, ny}, Cjx{nx, ny},
    Ey{nx, ny}, Jy{nx, ny}, Ceye{nx, ny}, Ceyh{nx, ny}, Cjy{nx, ny},
    Ez{nx, ny}, Jz{nx, ny}, Ceze{nx, ny}, Cezh{nx, ny}, Cjz{nx, ny},
    Hx{nx, ny}, Chxe{nx, ny}, Chxh{nx, ny},
    Hy{nx, ny}, Chye{nx, ny}, Chyh{nx, ny},
    Hz{nx, ny}, Chze{nx, ny}, Chzh{nx, ny}
  {
    DBG("EMData::EMData(nx, ny)", nx, ny);
  }
  
  ex_t Ex;
  ex_t Jx;
  ex_t Cexe;
  ex_t Cexh;
  ex_t Cjx;
  
  ey_t Ey;
  ey_t Jy;
  ey_t Ceye;
  ey_t Ceyh;
  ey_t Cjy;
  
  ez_t Ez;
  ez_t Jz;
  ez_t Ceze;
  ex_t Cezh;
  ez_t Cjz;
  
  hx_t Hx;
  hx_t Chxe;
  hx_t Chxh;
  
  hy_t Hy;
  hy_t Chye;
  hy_t Chyh;
  
  hz_t Hz;
  hz_t Chze;
  hz_t Chzh;
};




#endif //EM_DATA_H
