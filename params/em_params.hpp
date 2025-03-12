#ifndef EM_PARAMS_HPP
#define EM_PARAMS_HPP

using ex_t = double;
using ey_t = double;
using ez_t = double;
using hx_t = double;
using hy_t = double;
using hz_t = double;

// inline constexpr int BoundaryTypes[6] = {2, 2, 2, 2, 2, 2};


inline constexpr auto   PMLDepth = 10zu;
inline constexpr double PMLGrade = 3.0;
inline constexpr double PMLAlphaMax = 0.2;
inline constexpr double PMLKappaMax = 1.0;

#endif //EM_PARAMS_HPP
