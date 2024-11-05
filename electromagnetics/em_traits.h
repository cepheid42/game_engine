//
// Created by cepheid on 10/17/24.
//

#ifndef EM_TRAITS_H
#define EM_TRAITS_H

template<typename T>
concept FieldComponent = requires
{
  typename T::value_t;
  typename T::dimension_t;
  // typename T::array_t;
};


template<typename T, typename EMPTY>
concept is_empty_field = std::same_as<EMPTY, T>;

template<typename T, typename EMPTY>
concept is_1D_fields = requires
{
  requires is_empty_field<typename T::ex_t, EMPTY>;
  requires is_empty_field<typename T::ey_t, EMPTY>;
  requires is_empty_field<typename T::hx_t, EMPTY>;
  requires is_empty_field<typename T::hz_t, EMPTY>;

  requires !is_empty_field<typename T::ez_t, EMPTY>;
  requires !is_empty_field<typename T::hy_t, EMPTY>;
};

template<typename T, typename EMPTY>
concept is_TM_fields = requires
{
  requires is_empty_field<typename T::ex_t, EMPTY>;
  requires is_empty_field<typename T::ey_t, EMPTY>;
  requires is_empty_field<typename T::hz_t, EMPTY>;

  requires !is_empty_field<typename T::ez_t, EMPTY>;
  requires !is_empty_field<typename T::hx_t, EMPTY>;
  requires !is_empty_field<typename T::hy_t, EMPTY>;
};

#endif //EM_TRAITS_H
