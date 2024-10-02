//
// Created by cepheid on 9/23/24.
//

#ifndef EM_DATA_H
#define EM_DATA_H

#include "core/typelist.h"
#include "aydenstuff/array.h"

template<typename T, std::size_t... N>
struct EmptyArray {
  using value_t = T;
  using vector_t = std::vector<value_t>;
  using dimension_t = tf::tags::Dimension<sizeof...(N)>;

  EmptyArray() = default;
  explicit EmptyArray(std::size_t...) {}

  value_t operator()(std::size_t...) const { return static_cast<value_t>(0.0); }
};

template<typename T>
using EmptyArray1D = EmptyArray<T, 1>;

template<typename T>
using EmptyArray2D = EmptyArray<T, 2>;

template<typename T>
using EmptyArray3D = EmptyArray<T, 3>;


template<typename Array>
struct base_em_data {
  using array_t = Array;
  using value_t = typename Array::value_t;
  using dimension_t = typename Array::dimension_t;
};

template<typename Array>
struct ex_field {
  using ex_t = Array;
  // Array Ex;
  // Array Jx;
  // Array Cexe;
  // Array Cexh;
};

template<typename Array>
struct ey_field {
  Array Ey;
  Array Jy;
  Array Ceye;
  Array Ceyh;
};

template<typename Array>
struct ez_field {
  Array Ez;
  Array Jz;
  Array Ceze;
  Array Cezh;
};

template<typename Base, int D>
struct Discriminator : Base {};

template<typename ExSetter, typename EySetter, typename EzSetter>
struct FieldSelector
: Discriminator<ExSetter, 1>, Discriminator<EySetter, 2>, Discriminator<EzSetter, 3>
{};

// template<typename Array>
struct DefaultFields {
  using ex_t = base_em_data<EmptyArray1D<double>>;
  using ey_t = base_em_data<EmptyArray1D<double>>;
  using ez_t = base_em_data<EmptyArray1D<double>>;
};

// template<typename Array>
struct DefaultFieldArgs : virtual DefaultFields {};

template<typename Field>
struct exfield_is : virtual DefaultFields {
  using ex_t = Field;
};

template<typename Field>
struct eyfield_is : virtual DefaultFields {
  using ey_t = Field;
};
template<typename Field>
struct ezfield_is : virtual DefaultFields {
  using ez_t = Field;
};

template<typename Fields>
struct EMData {
  using ex_t = typename Fields::Types::ex_t;
  using ey_t = typename Fields::Types::ey_t;
  using ez_t = typename Fields::Types::ez_t;

  ex_t Ez;
  ex_t Jz;
  ex_t Ceze;
  ex_t Cezh;
};

struct emdataNone {
  using array_t = tf::types::Array1D<double>;
  using Types = FieldSelector<
    DefaultFieldArgs,
    DefaultFieldArgs,
    DefaultFieldArgs
  >;
};

struct emdata1D {
  using array_t = tf::types::Array1D<double>;
  using Types = FieldSelector<
    DefaultFieldArgs,
    DefaultFieldArgs,
    ezfield_is<ez_field<array_t>>
  >;
};

using FieldTypeList = TypeList<emdataNone, emdata1D>;

using emdata_t = EMData<TypeListAt<1, FieldTypeList>>;

// template<typename T>
// struct em_data_1d : base_em_data<tf::types::Array1D<T>> {
//   using array_t = typename base_em_data<tf::types::Array1D<T>>::array_t;
//   using dimension_t = typename base_em_data<tf::types::Array1D<T>>::dimension_t;
//   using value_t = typename base_em_data<tf::types::Array1D<T>>::value_t;
//
//   em_data_1d() = default;
//
//   em_data_1d(const size_t nx)
//   : Ez{nx},
//     Hy{nx - 1},
//     Jz{nx},
//     Ceze{nx},
//     Cezh{nx},
//     Chye{nx - 1},
//     Chyh{nx - 1},
//     Cjz{nx}
//   {}
//
//   array_t Ez;
//   array_t Hy;
//   array_t Jz;
//   array_t Ceze;
//   array_t Cezh;
//   array_t Chye;
//   array_t Chyh;
//   array_t Cjz;
// };
//
// template<typename Array>
// struct em_data_tm : base_em_data<Array> {
//   using array_t = typename base_em_data<Array>::array_t;
//   using dimension_t = typename base_em_data<Array>::dimension_t;
//   using value_t = typename base_em_data<Array>::value_t;
//   using empty_t = EmptyArray2D<value_t>;
//
//   em_data_tm() = default;
//
//   em_data_tm(const size_t nx, const size_t ny)
//   : base_em_data<Array>(),
//     Ez{nx, ny},
//     Hx{nx, ny - 1},
//     Hy{nx - 1, ny},
//     Jz{nx, ny},
//     Ceze{nx, ny},
//     Cezh{nx, ny},
//     Chxe{nx, ny - 1},
//     Chxh{nx, ny - 1},
//     Chye{nx - 1, ny},
//     Chyh{nx - 1, ny},
//     Cjz{nx, ny}
//   {}
//   empty_t Ex{};
//   empty_t Ey{};
//   array_t Ez;
//   array_t Hx;
//   array_t Hy;
//   empty_t Hz{};
//   empty_t Jx{};
//   empty_t Jy{};
//   array_t Jz;
//   empty_t Cexe{};
//   empty_t Cexh{};
//   empty_t Ceye{};
//   empty_t Ceyh{};
//   array_t Ceze;
//   array_t Cezh;
//   array_t Chxe;
//   array_t Chxh;
//   array_t Chye;
//   array_t Chyh;
//   empty_t Chze{};
//   empty_t Chzh{};
//   empty_t Cjx{};
//   empty_t Cjy{};
//   array_t Cjz;
// };
//
// template<typename Array>
// struct em_data_te : base_em_data<Array> {
//   using array_t = typename base_em_data<Array>::array_t;
//   using dimension_t = typename base_em_data<Array>::dimension_t;
//   using value_t = typename base_em_data<Array>::value_t;
//
//   em_data_te() = default;
//
//   em_data_te(const size_t nx, const size_t ny)
//   : base_em_data<Array>(),
//     Ex{nx - 1, ny},
//     Ey{nx, ny - 1},
//     Hz{nx - 1, ny - 1},
//     Jx{nx - 1, ny},
//     Jy{nx, ny - 1},
//     Cexe{nx - 1, ny},
//     Cexh{nx - 1, ny},
//     Ceye{nx, ny - 1},
//     Ceyh{nx, ny - 1},
//     Chze{nx - 1, ny - 1},
//     Chzh{nx - 1, ny - 1},
//     Cjx{nx - 1, ny},
//     Cjy{nx, ny - 1}
//   {}
//   array_t Ex;
//   array_t Ey;
//   array_t Hz;
//   array_t Jx;
//   array_t Jy;
//   array_t Cexe;
//   array_t Cexh;
//   array_t Ceye;
//   array_t Ceyh;
//   array_t Chze;
//   array_t Chzh;
//   array_t Cjx;
//   array_t Cjy;
// };
//
// template<typename Array>
// struct em_data_3d : base_em_data<Array> {
//   using array_t = typename base_em_data<Array>::array_t;
//   using dimension_t = typename base_em_data<Array>::dimension_t;
//   using value_t = typename base_em_data<Array>::value_t;
//
//   em_data_3d() = default;
//
//   em_data_3d(const size_t nx, const size_t ny, const size_t nz)
//   : base_em_data<Array>(),
//     Ex{nx - 1, ny, nz},
//     Ey{nx, ny - 1, nz},
//     Ez{nx, ny, nz - 1},
//     Hx{nx, ny - 1, nz - 1},
//     Hy{nx - 1, ny, nz - 1},
//     Hz{nx - 1, ny - 1, nz},
//     Jx{nx - 1, ny, nz},
//     Jy{nx, ny - 1, nz},
//     Jz{nx, ny, nz - 1},
//     Cexe{nx - 1, ny, nz},
//     Cexh{nx - 1, ny, nz},
//     Ceye{nx, ny - 1, nz - 1},
//     Ceyh{nx, ny - 1, nz - 1},
//     Ceze{nx, ny, nz - 1},
//     Cezh{nx, ny, nz - 1},
//     Chxe{nx, ny - 1, nz - 1},
//     Chxh{nx, ny - 1, nz - 1},
//     Chye{nx - 1, ny, nz - 1},
//     Chyh{nx - 1, ny, nz - 1},
//     Chze{nx - 1, ny - 1, nz},
//     Chzh{nx - 1, ny - 1, nz},
//     Cjx{nx - 1, ny, nz},
//     Cjy{nx, ny - 1, nz},
//     Cjz{nx, ny, nz - 1}
//   {}
//
//   array_t Ex;
//   array_t Ey;
//   array_t Ez;
//   array_t Hx;
//   array_t Hy;
//   array_t Hz;
//   array_t Jx;
//   array_t Jy;
//   array_t Jz;
//   array_t Cexe;
//   array_t Cexh;
//   array_t Ceye;
//   array_t Ceyh;
//   array_t Ceze;
//   array_t Cezh;
//   array_t Chxe;
//   array_t Chxh;
//   array_t Chye;
//   array_t Chyh;
//   array_t Chze;
//   array_t Chzh;
//   array_t Cjx;
//   array_t Cjy;
//   array_t Cjz;
// };

#endif //EM_DATA_H
