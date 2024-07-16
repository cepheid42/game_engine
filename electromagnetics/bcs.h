//
// Created by cepheid on 7/16/24.
//

#ifndef BCS_H
#define BCS_H

#include "enums.h"
#include "curls.h"
#include "array.h"

struct PMLBoundary2D
{
    Array2D psi;
    Array2D c_psi;
    Array2D c_f;
};

// struct BCS
// {
//     PMLBoundary2D ez_pml;
//     PMLBoundary2D hx_pml;
//     PMLBoundary2D hy_pml;
// };

template<typename TL, FieldType FT, Derivative CURL, typename... IDXS>
struct PMLConvolutionFunctor {
    using curl1 = curl<CURL, FT, IDXS...>;
    using F = TypeListAt<0, TL>;
    using P = TypeListAt<1, TL>;
    using D1 = TypeListAt<2, TL>;
    using CF = TypeListAt<3, TL>;
    using CP = TypeListAt<4, TL>;
    using CD = TypeListAt<5, TL>;

    static auto apply(F& psi, const D1& d, const CP& c_psi, const CD& c_d, IDXS... idxs) {
        const auto prev = c_psi(idxs...) * psi(idxs...);
        const auto diff = c_d(idxs...) * curl1::apply(d, idxs...);
        psi(idxs...) = prev + diff;
    }
};


template<typename TL, FieldType FT, Derivative CURL>
struct PMLConvIntegrator2D {
    using F = TypeListAt<0, TL>;
    using P = TypeListAt<1, TL>;
    using D1 = TypeListAt<2, TL>;
    using CF = TypeListAt<3, TL>;
    using CP = TypeListAt<4, TL>;
    using CD = TypeListAt<5, TL>;
    using conv_func = PMLConvolutionFunctor<TL, FT, CURL, size_t, size_t>;

    static auto apply(P& psi, const D1& d, const CP& c_psi, const CD& c_d, const IntegratorOffsets<2>& o) {
        const auto [x0, x1, y0, y1] = o.offsets;
        for (size_t i = x0; i < psi.shape[0] - x1; i++) {
            for (size_t j = y0; j < psi.shape[1] - y1; j++) {
                conv_func::apply(psi, d, c_psi, c_d, i, j);
            }
        }
    }
};

template<typename TL, FieldType FT, typename... IDXS>
struct PMLUpdateFunctor {
    using F = TypeListAt<0, TL>;
    using P = TypeListAt<1, TL>;
    using D1 = TypeListAt<2, TL>;
    using CF = TypeListAt<3, TL>;
    using CP = TypeListAt<4, TL>;
    using CD = TypeListAt<5, TL>;

    static auto apply(F& f, const P& psi, const CF& c_f, const IDXS... idxs) {

    }
};

template<typename TL, FieldType FT>
struct PMLUpdateIntegrator2D {
    using F = TypeListAt<0, TL>;
    using P = TypeListAt<1, TL>;
    using D1 = TypeListAt<2, TL>;
    using CF = TypeListAt<3, TL>;
    using CP = TypeListAt<4, TL>;
    using CD = TypeListAt<5, TL>;
    using update_func = PMLUpdateFunctor<TL, FT,  size_t, size_t>;

    static auto apply(F& f, const P& psi, const CF& c_f, const IntegratorOffsets<2>& o) {
        const auto [x0, x1, y0, y1] = o.offsets;
        for (size_t i = x0; i < f.shape[0] - x1; i++) {
            for (size_t j = y0; j < f.shape[1] - y1; j++) {
                update_func::apply(f, psi, c_f, i, j);
            }
        }
    }
};

// template<>
// struct PMLBoundary
// {
//     static void updateE(auto& pml, auto& field, const auto& coeff) {
//         CI::apply(pml.psi, field, coeff);
//     }
//
//     static void updateB() {}
// };


#endif //BCS_H
