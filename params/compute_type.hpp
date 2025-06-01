#ifndef COMPUTE_TYPE_HPP
#define COMPUTE_TYPE_HPP

using compute_t = double;

constexpr compute_t operator""_fp(const long double x) {
   return static_cast<compute_t>(x);
}

#endif //COMPUTE_TYPE_HPP
