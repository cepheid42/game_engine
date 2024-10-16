//
// Created by cepheid on 9/23/24.
//

#ifndef EM_DATA_H
#define EM_DATA_H

#include "core/typelist.h"
#include "aydenstuff/array.h"
#include "core/debug.h"

using tf::types::Array1D;
using tf::types::Array2D;

template<typename T>
concept FieldComponent = requires
{
  typename T::arr_t;
};



template<typename T, std::size_t N>
struct EmptyArray {
  using value_t = T;
  using vector_t = std::vector<value_t>;
  using dimension_t = tf::tags::Dimension<N>;

  EmptyArray() = default;
  explicit EmptyArray(std::size_t...) {}

  constexpr value_t operator[](std::size_t) const { return static_cast<value_t>(0.0); }
  constexpr value_t operator()(std::size_t...) const { return static_cast<value_t>(0.0); }
};

template<typename T>
using EmptyArray1D = EmptyArray<T, 1>;

template<typename T>
using EmptyArray2D = EmptyArray<T, 2>;

template<typename T>
using EmptyArray3D = EmptyArray<T, 3>;

template<typename T>
concept is_empty_field = std::same_as<EmptyArray<typename T::value_t, T::dimension_t::value>, T>;

template<typename T>
concept is_1D_fields = requires
{
  requires is_empty_field<typename T::ex_t>;
  requires is_empty_field<typename T::ey_t>;
  requires is_empty_field<typename T::hx_t>;
  requires is_empty_field<typename T::hz_t>;
  
  requires !is_empty_field<typename T::ez_t>;
  requires !is_empty_field<typename T::hy_t>;
};

template<typename T>
concept is_TM_fields = requires
{
  requires is_empty_field<typename T::ex_t>;
  requires is_empty_field<typename T::ey_t>;
  requires is_empty_field<typename T::hz_t>;

  requires !is_empty_field<typename T::ez_t>;
  requires !is_empty_field<typename T::hx_t>;
  requires !is_empty_field<typename T::hy_t>;
};

template<typename Array>
struct enable_field {
  using value_t = typename Array::value_t;
  using dimension_t = typename Array::dimension_t;
  using arr_t = Array;
};

template<typename Array>
struct disable_field {
  using value_t = typename Array::value_t;
  using dimension_t = typename Array::dimension_t;
  using arr_t = EmptyArray<value_t, dimension_t::value>;
};

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


template<typename T>
using emdataNone = EMData<
  disable_field<EmptyArray1D<T>>, // Ex
  disable_field<EmptyArray1D<T>>, // Ey
  disable_field<EmptyArray1D<T>>, // Ez
  disable_field<EmptyArray1D<T>>, // Hx
  disable_field<EmptyArray1D<T>>, // Hy
  disable_field<EmptyArray1D<T>>  // Hz
>;

template<typename T>
using emdata1D = EMData<
  disable_field<EmptyArray1D<T>>, // Ex
  disable_field<EmptyArray1D<T>>, // Ey
  enable_field<Array1D<T>>,       // Ez
  disable_field<Array1D<T>>,      // Hx
  enable_field<Array1D<T>>,       // Hy
  disable_field<Array1D<T>>       // Hz
>;

template<typename T>
using emdataTM = EMData<
  disable_field<EmptyArray2D<T>>, // Ex
  disable_field<EmptyArray2D<T>>, // Ey
  enable_field<Array2D<T>>,       // Ez
  enable_field<Array2D<T>>,       // Hx
  enable_field<Array2D<T>>,       // Hy
  disable_field<Array2D<T>>       // Hz
>;

template<typename T>
using EMDataTL = TypeList<emdataNone<T>, emdata1D<T>, emdataTM<T>>;

template<typename T>
using emdata_t = TypeListAt<2, EMDataTL<T>>;

#endif //EM_DATA_H
