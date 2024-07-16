//
// Created by cepheid on 7/16/24.
//

#ifndef ENUMS_H
#define ENUMS_H

enum class Derivative{ DX, DY, DZ, NoOp };
enum class FieldType{ E, H };
enum class FieldComp{ X, Y, Z };

template<size_t DIM>
struct IntegratorOffsets {
    std::array<size_t, 2 * DIM> offsets;
};

#endif //ENUMS_H
