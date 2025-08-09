#ifndef EM_PARAMS_HPP
#define EM_PARAMS_HPP

#include <array>

enum class EMFace { X, Y, Z };
enum class EMSide { Lo, Hi };

// Reflecting = 0, Periodic = 1, PML = 2
inline constexpr std::array<std::size_t, 6> BCSelect = {0, 0, 0, 0, 0, 0};

#endif //EM_PARAMS_HPP
