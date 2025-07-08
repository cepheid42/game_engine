#ifndef MORTON_HPP
#define MORTON_HPP

#include <array>

namespace tf {
// HELPER METHOD: Magic bits encoding (helper method)
static constexpr std::size_t morton3D_SplitBy3bits(std::size_t x) {
   x = x & 0x1fffff;
   x = (x | x << 32) & 0x1f00000000ffff;
   x = (x | x << 16) & 0x1f0000ff0000ff;
   x = (x | x << 8) & 0x100f00f00f00f00f;
   x = (x | x << 4) & 0x10c30c30c30c30c3;
   x = (x | x << 2) & 0x1249249249249249;
   return x;
}

// ENCODE 3D Morton code : Magic bits method
// This method uses certain bit patterns (magic bits) to split bits in the coordinates
constexpr std::size_t morton_encode(const std::size_t x, const std::size_t y, const std::size_t z) {
   return morton3D_SplitBy3bits(z) | (morton3D_SplitBy3bits(y) << 1) | (morton3D_SplitBy3bits(x) << 2);
}

constexpr std::size_t morton_encode(const auto& vec) {
   return morton3D_SplitBy3bits(vec[2]) | (morton3D_SplitBy3bits(vec[1]) << 1) | (morton3D_SplitBy3bits(vec[0]) << 2);
}


// HELPER METHOD for Magic bits decoding
static constexpr std::size_t morton3D_GetThirdBits(const std::size_t m) {
   std::size_t x = m & 0x1249249249249249;
   x = (x ^ (x >> 2)) & 0x10c30c30c30c30c3;
   x = (x ^ (x >> 4)) & 0x100f00f00f00f00f;
   x = (x ^ (x >> 8)) & 0x1f0000ff0000ff;
   x = (x ^ (x >> 16)) & 0x1f00000000ffff;
   x = (x ^ (x >> 32)) & 0x1fffff;
   return x;
}

// DECODE 3D Morton code : Magic bits
// This method splits the morton codes bits by using certain patterns (magic bits)
constexpr std::array<std::size_t, 3> morton_decode(const std::size_t m) {
   return {
      morton3D_GetThirdBits(m >> 2),
      morton3D_GetThirdBits(m >> 1),
      morton3D_GetThirdBits(m)
   };
}
} // end namespace tf
#endif //MORTON_HPP
