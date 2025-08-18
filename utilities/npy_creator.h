#ifndef GAME_ENGINE_NPY_CREATOR_H
#define GAME_ENGINE_NPY_CREATOR_H

#include "array.hpp"
#include "vec3.hpp"

#include <string>
#include <print>

using namespace std::string_literals;

template<bool> consteval auto getTypeStr() { return "b"s; }
template<typename> consteval auto getTypeStr() { return "?"s; }
template<std::floating_point> consteval auto getTypeStr() { return "f"s; }
template<std::signed_integral> consteval auto getTypeStr() { return "i"s; }
template<std::unsigned_integral> consteval auto getTypeStr() { return "u"s; }

template<typename T>
void create_npy_header(const auto& shape) {
   const auto format = "0x93NUMPY0x020x00{}{}\n"s;

   auto dict = "{'descr': "s +
               "'<"s + getTypeStr<T>() + std::to_string(sizeof(T)) + "', "s +
               "'fortran_order': False, "s +
               "'shape': ("s + std::to_string(shape[0]) + ", "s + std::to_string(shape[1]) + ", "s + std::to_string(shape[2]) + ")}"s;

   std::print("HEADER_LEN = {}", dict.size());
}


#endif //GAME_ENGINE_NPY_CREATOR_H