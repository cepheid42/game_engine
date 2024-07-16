//
// Created by cepheid on 7/15/24.
//

#ifndef REGULAR_EM_H
#define REGULAR_EM_H

using fptype = double;

template<typename F, typename D1, typename D2, typename C1, typename C2>
void update_Ez2D(F& Ez, const D1& Hy, const D2& Hx, const C1& c_f, const C2& c_d) {
    for (size_t i = 1; i < Ez.shape[0] - 1; i++) {
        for (size_t j = 1; j < Ez.shape[1] - 1; j++) {
            const auto self = c_f(i, j) * Ez(i, j);
            const auto diff1 = Hy(i, j) - Hy(i - 1, j);
            const auto diff2 = Hx(i, j) - Hx(i, j - 1);
            const auto diff = c_d(i, j) * (diff1 - diff2);
            Ez(i, j) = self + diff;
        }
    }
}

template<typename F, typename D1, typename C1, typename C2>
void update_Hx2D(F& Hx, const D1& Ez, const C1& c_f, const C2& c_d) {
    for (size_t i = 0; i < Hx.shape[0]; i++) {
        for (size_t j = 0; j < Hx.shape[1]; j++) {
            const auto self = c_f(i, j) * Hx(i, j);
            const auto diff1 = Ez(i, j + 1) - Ez(i, j);
            const auto diff = c_d(i, j) * diff1;
            Hx(i, j) = self - diff;
        }
    }
}

template<typename F, typename D1, typename C1, typename C2>
void update_Hy2D(F& Hy, const D1& Ez, const C1& c_f, const C2& c_d) {
    for (size_t i = 0; i < Hy.shape[0]; i++) {
        for (size_t j = 0; j < Hy.shape[1]; j++) {
            const auto self = c_f(i, j) * Hy(i, j);
            const auto diff1 = Ez(i + 1, j) - Ez(i, j);
            const auto diff = c_d(i, j) * diff1;
            Hy(i, j) = self + diff;
        }
    }
}

template<typename F, typename D1, typename D2, typename C1, typename C2>
void update_Ex3D(F& Ex, const D1& Hz, const D2& Hy, const C1& c_f, const C2& c_d) {
    for (size_t i = 0; i < Ex.shape[0]; i++) {
        for (size_t j = 1; j < Ex.shape[1] - 1; j++) {
            for (size_t k = 1; k < Ex.shape[2] - 1; k++) {
                const auto self = c_f(i, j, k) * Ex(i, j, k);
                const auto diff1 = Hz.backward_diff_y(i, j, k); // Hz(i, j, k) - Hz(i, j - 1, k);
                const auto diff2 = Hy.backward_diff_z(i, j, k); // Hy(i, j, k) - Hy(i, j, k - 1);
                const auto diff = c_d(i, j, k) * (diff1 - diff2);
                Ex(i, j, k) = self + diff;
            }
        }
    }
}

template<typename F, typename D1, typename D2, typename C1, typename C2>
void update_Ey3D(F& Ey, const D1& Hx, const D2& Hz, const C1& c_f, const C2& c_d) {
    for (size_t i = 1; i < Ey.shape[0] - 1; i++) {
        for (size_t j = 0; j < Ey.shape[1]; j++) {
            for (size_t k = 1; k < Ey.shape[2] - 1; k++) {
                const auto self = c_f(i, j, k) * Ey(i, j, k);
                const auto diff1 = Hx.backward_diff_z(i, j, k); //Hx(i, j, k) - Hx(i, j, k - 1);
                const auto diff2 = Hz.backward_diff_x(i, j, k);//Hz(i, j, k) - Hz(i - 1, j, k);
                const auto diff = c_d(i, j, k) * (diff1 - diff2);
                Ey(i, j, k) = self + diff;
            }
        }
    }
}

template<typename F, typename D1, typename D2, typename C1, typename C2>
void update_Ez3D(F& Ez, const D1& Hy, const D2& Hx, const C1& c_f, const C2& c_d) {
    for (size_t i = 1; i < Ez.shape[0] - 1; i++) {
        for (size_t j = 1; j < Ez.shape[1] - 1; j++) {
            for (size_t k = 0; k < Ez.shape[2]; k++) {
                const auto self = c_f(i, j, k) * Ez(i, j, k);
                const auto diff1 = Hy.backward_diff_x(i, j, k); //Hy(i, j, k) - Hy(i - 1, j, k);
                const auto diff2 = Hx.backward_diff_y(i, j, k); //Hx(i, j, k) - Hx(i, j - 1, k);
                const auto diff = c_d(i, j, k) * (diff1 - diff2);
                Ez(i, j, k) = self + diff;
            }
        }
    }
}

template<typename F, typename D1, typename D2, typename C1, typename C2>
void update_Hx3D(F& Hx, const D1& Ey, const D2& Ez, const C1& c_f, const C2& c_d) {
    for (size_t i = 0; i < Hx.shape[0]; i++) {
        for (size_t j = 0; j < Hx.shape[1]; j++) {
            for (size_t k = 0; k < Hx.shape[2]; k++) {
                const auto self = c_f(i, j, k) * Hx(i, j, k);
                const auto diff1 = Ey.forward_diff_z(i, j, k); //Ey(i, j, k + 1) - Ey(i, j, k);
                const auto diff2 = Ez.forward_diff_y(i, j, k); //Ez(i, j + 1, k) - Ez(i, j, k);
                const auto diff = c_d(i, j, k) * (diff1 - diff2);
                Hx(i, j, k) = self + diff;
            }
        }
    }
}

template<typename F, typename D1, typename D2, typename C1, typename C2>
void update_Hy3D(F& Hy, const D1& Ez, const D2& Ex, const C1& c_f, const C2& c_d) {
    for (size_t i = 0; i < Hy.shape[0]; i++) {
        for (size_t j = 0; j < Hy.shape[1]; j++) {
            for (size_t k = 0; k < Hy.shape[2]; k++) {
                const auto self = c_f(i, j, k) * Hy(i, j, k);
                const auto diff1 = Ez.forward_diff_x(i, j, k); //Ez(i + 1, j, k) - Ez(i, j, k);
                const auto diff2 = Ex.forward_diff_z(i, j, k); //Ex(i, j, k + 1) - Ex(i, j, k);
                const auto diff = c_d(i, j, k) * (diff1 - diff2);
                Hy(i, j, k) = self + diff;
            }
        }
    }
}

template<typename F, typename D1, typename D2, typename C1, typename C2>
void update_Hz3D(F& Hz, const D1& Ex, const D2& Ey, const C1& c_f, const C2& c_d) {
    for (size_t i = 0; i < Hz.shape[0]; i++) {
        for (size_t j = 0; j < Hz.shape[1]; j++) {
            for (size_t k = 0; k < Hz.shape[2]; k++) {
                const auto self = c_f(i, j, k) * Hz(i, j, k);
                const auto diff1 = Ex.forward_diff_y(i, j, k); //Ex(i, j + 1, k) - Ex(i, j, k);
                const auto diff2 = Ey.forward_diff_x(i, j, k); //Ey(i + 1, j, k) - Ey(i, j, k);
                const auto diff = c_d(i, j, k) * (diff1 - diff2);
                Hz(i, j, k) = self + diff;
            }
        }
    }
}


#endif //REGULAR_EM_H
