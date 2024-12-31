//
// Created by cepheid on 10/17/24.
//

#ifndef EM_TRAITS_H
#define EM_TRAITS_H

namespace tf::electromagnetics::traits
{
  enum class EMFace { X, Y, Z };
  enum class EMSide { Lo, Hi };
  enum class EMComponent { E, H };
  enum class Derivative { DX, DY, DZ, NoOp };

  template<EMFace F, EMSide S>
  struct periodic_t {
    static constexpr EMFace face = F;
    static constexpr EMSide side = S;
  };

  template<EMFace F, EMSide S>
  struct pml_t {
    static constexpr EMFace face = F;
    static constexpr EMSide side = S;
  };

  template<typename T>
  concept FieldComponent = requires (T t)
  {
    typename T::value_t;
    typename T::dimension_t;

    { t.nx() }-> std::same_as<size_t>;
    { t.ny() }-> std::same_as<size_t>;
    { t.nz() }-> std::same_as<size_t>;
  };

  // template<typename T, typename EMPTY>
  // concept is_empty_field = std::same_as<T, EMPTY>;

  template <template <typename...> typename C, typename   T> struct type_instance_impl : std::false_type {};
  template <template <typename...> typename C, typename...T> struct type_instance_impl<C, C<T...>> : std::true_type {};
  template <template <typename...> typename C, typename   T> concept instance_of_type = type_instance_impl<C, std::decay_t<T>>::value;


  template<typename T>
  concept is_periodic = std::derived_from<T, periodic_t<T::face, T::side>>;

  template<typename T>
  concept is_pml = std::derived_from<T, pml_t<T::face, T::side>>;
  //
} // end namespace tf::electromagnetics::traits
#endif //EM_TRAITS_H
