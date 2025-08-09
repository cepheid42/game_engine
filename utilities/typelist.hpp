#ifndef TYPELIST_HPP
#define TYPELIST_HPP

#include <cstddef>

//===== TypeList =====
template<typename...>
struct TypeList{};

//===== TypeListSize =====
template<typename> struct TypeListSizeImpl; // declares template for implementation

template<typename... Types>
struct TypeListSizeImpl<TypeList<Types...>> {
   // implementation for generic TypeList
   static constexpr std::size_t value = sizeof...(Types);
};

template<typename Types>
constexpr std::size_t TypeListSize = TypeListSizeImpl<Types>::value; // implementation tied to interface

//===== TypeListAt =====
template<std::size_t, typename> struct TypeListAtImpl; // declares template for instatiation

template<std::size_t I, typename Type, typename... Types>
struct TypeListAtImpl<I, TypeList<Type, Types...>> {
   // implementation for generic Typelist
   using type = TypeListAtImpl<I - 1, TypeList<Types...>>::type; // peels off first element of TypeList
};

template<typename Type, typename... Types>
struct TypeListAtImpl<0, TypeList<Type, Types...>> {
   // implementation for I = 0
   using type = Type;
};

template<std::size_t I, typename Types>
using TypeListAt = TypeListAtImpl<I, Types>::type; // implementation tied to interface

#endif //TYPELIST_HPP