#ifndef TRAITS_HPP
#define TRAITS_HPP

#include <concepts>

struct periodic_t {};
struct pml_t {};
struct null_t {};

template<typename T> concept is_periodic = std::derived_from<T, periodic_t>;
template<typename T> concept is_pml = std::derived_from<T, pml_t>;
template<typename T> concept is_null = std::derived_from<T, null_t>;

#endif //TRAITS_HPP
