//
// Created by cepheid on 7/16/24.
//

#ifndef BCS_H
#define BCS_H

#include "enums.h"
#include "curls.h"
#include "array.h"

//=================== PML Functors ========================
//=========================================================
template<size_t Axis, size_t Side>
size_t getPMLIndexOffset(const size_t dim[2], const size_t d_pml) {
    if constexpr (Side == 0) {
        return 0; // X0/Y0/Z0
    }

    if constexpr(Axis == 0) {
        return dim[1] - d_pml; // X1
    } else {
        return dim[1] * (dim[0] - d_pml); // Y1
    }
}

template<size_t DIM>
struct PMLOffsets {
    std::array<size_t, 2 * DIM + 1> offsets;
};


template<typename TL, FieldType FT, Derivative CURL, typename... IDXS>
struct PMLConvolutionFunctor {
    using curl1 = curl<CURL, FT, IDXS...>;
    using F = TypeListAt<0, TL>;
    using P = TypeListAt<1, TL>;
    using D1 = TypeListAt<2, TL>;
    using CF = TypeListAt<3, TL>;
    using CP = TypeListAt<4, TL>;
    using CD = TypeListAt<5, TL>;

    static auto apply(P& psi, const D1& d, const CP& c_psi, const CD& c_d, const size_t f_offset, const IDXS... idxs) {
        const auto prev = c_psi(idxs...) * psi(idxs...);
        const auto diff = c_d(idxs..., f_offset) * curl1::apply(d, idxs..., f_offset);
        psi(idxs...) = prev + diff;
    }
};

template<typename TL, FieldComp FC, typename... IDXS>
struct PMLUpdateFunctor {
    using F = TypeListAt<0, TL>;
    using P = TypeListAt<1, TL>;
    using D1 = TypeListAt<2, TL>;
    using CF = TypeListAt<3, TL>;
    using CP = TypeListAt<4, TL>;
    using CD = TypeListAt<5, TL>;

    static auto apply(F& f, const P& psi, const CF& c_f, const size_t f_offset, const IDXS... idxs) {
        if constexpr (FC == FieldComp::X) {
            f(idxs..., f_offset) += c_f(idxs..., f_offset) * psi(idxs...);
        } else {
            f(idxs..., f_offset) -= c_f(idxs..., f_offset) * psi(idxs...);
        }
    }
};


template<typename TL, FieldType FT, FieldComp FC, Derivative CURL>
struct PMLIntegrator2D {
    using F = TypeListAt<0, TL>;
    using P = TypeListAt<1, TL>;
    using D1 = TypeListAt<2, TL>;
    using CF = TypeListAt<3, TL>;
    using CP = TypeListAt<4, TL>;
    using CD = TypeListAt<5, TL>;
    using conv_func = PMLConvolutionFunctor<TL, FT, CURL, size_t, size_t>;
    using update_func = PMLUpdateFunctor<TL, FC, size_t, size_t>;

    static auto apply(F& f, P& psi, const D1& d, const CF& c_f, const CP& c_psi, const CD& c_d, const PMLOffsets<2>& o) {
        const auto [x0, x1, y0, y1, lo] = o.offsets;

        for (size_t i = x0; i < psi.shape[0] - x1; i++) {
            for (size_t j = y0; j < psi.shape[1] - y1; j++) {
                conv_func::apply(psi, d, c_psi, c_d, lo, i, j);
                update_func::apply(f, psi, c_f, lo, i, j);
            }
        }
    }
};

//=================== Boundary Functors ========================
//==============================================================
struct NullBoundary {
    static constexpr auto apply(...) {}
};

template<typename TL>
struct PeriodicBoundary {
    static auto apply(...) {}
};

//=================== Boundary Condition Class ========================
//=====================================================================
template<typename EBC, typename HBC>
struct EMBoundary {
    static void updateE(auto& f, auto& bcs, const auto& c) { EBC::apply(f.Ex, f.Ey, f.Ez); }
    static void updateH() { HBC::apply(); }
};




#endif //BCS_H
