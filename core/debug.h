//
// Created by akis on 1/29/24.
// References:
// https://devblogs.microsoft.com/cppblog/build-throughput-series-more-efficient-template-metaprogramming/
// https://www.codingwiththomas.com/blog/getting-started-with-typelists
//

#include <cstddef>

#ifndef TRIFORCE_TYPELIST_H
#define TRIFORCE_TYPELIST_H

using std::size_t;

//===== TypeList =====
template<typename... Types>
struct TypeList{};

//===== TypeListSize =====
template<typename> struct TypeListSizeImpl; // declares template for implementation

template<typename... Types>
struct TypeListSizeImpl<TypeList<Types...>> // implementation for generic TypeList
{
  static constexpr size_t value = sizeof...(Types);
};

template<typename Types>
constexpr size_t TypeListSize = TypeListSizeImpl<Types>::value; // implementation tied to interface

//===== TypeListAt =====
template<size_t, typename> struct TypeListAtImpl; // declares template for instatiation

template<size_t I, typename Type, typename... Types>
struct TypeListAtImpl<I, TypeList<Type, Types...>> // implementation for generic Typelist
{
  using type = typename TypeListAtImpl<I - 1, TypeList<Types...>>::type; // peels off first element of TypeList
};

template<typename Type, typename... Types>
struct TypeListAtImpl<0, TypeList<Type, Types...>> // implementation for I = 0
{
  using type = Type;
};

template<size_t I, typename Types>
using TypeListAt = typename TypeListAtImpl<I, Types>::type; // implementation tied to interface

// //===== Conditional Type Check =====
// template<bool COND, typename TrueType, typename FalseType>
// struct IfThenElseImpl {
//   using type = TrueType;
// };
//
// template<typename TrueType, typename FalseType>
// struct IfThenElseImpl<false, TrueType, FalseType> {
//   using type = FalseType;
// };
//
// template<bool COND, typename TrueType, typename FalseType>
// using IfThenElse = typename IfThenElseImpl<COND, TrueType, FalseType>::type;



// //===== Scattered Heirarchy =====
// class EmptyType {};
//
// template <typename TList, template <typename> typename Unit>
// class GenScatterHierarchy;
//
// template <typename Head, typename Tail, template <typename> typename Unit>
// class GenScatterHierarchy<TypeList<Head, Tail>, Unit>
// : public GenScatterHierarchy<Head, Unit>
// , public GenScatterHierarchy<Tail, Unit>
// {
// public:
//   typedef typename TypeList<Head, Tail> TList;
//   typedef typename GenScatterHierarchy<Head, Unit> LeftBase;
//   typedef typename GenScatterHierarchy<Tail, Unit> RightBase;
// };
// // Pass an atomic type (non-typelist) to Unit
// template <typename AtomicType, template <typename> typename Unit>
// class GenScatterHierarchy : public Unit<AtomicType>
// {
//   typedef typename Unit<AtomicType> LeftBase;
// };
// // Do nothing for NullType
// template <template <typename> typename Unit>
// class GenScatterHierarchy<EmptyType, Unit> {};


#endif  // TRIFORCE_TYPELIST_H
