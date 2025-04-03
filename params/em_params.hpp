#ifndef EM_PARAMS_HPP
#define EM_PARAMS_HPP

enum class BCType { PML, Periodic };

inline constexpr auto BoundaryType = BCType::Periodic;

inline constexpr auto PMLDepth = 10zu;
inline constexpr auto PMLGrade = 3.0;
inline constexpr auto PMLAlphaMax = 0.2;
inline constexpr auto PMLKappaMax = 1.0;

#endif //EM_PARAMS_HPP
