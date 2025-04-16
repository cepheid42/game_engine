#ifndef EM_PARAMS_HPP
#define EM_PARAMS_HPP

enum class BCType { PML, Periodic };

inline constexpr auto BoundaryType = BCType::PML;

inline constexpr auto PMLDepth = 50zu;
inline constexpr auto PMLGrade = 3.5;
inline constexpr auto PMLAlphaMax = 0.2;
inline constexpr auto PMLKappaMax = 1.0;

#endif //EM_PARAMS_HPP
