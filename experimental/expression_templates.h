//
// Created by cepheid on 8/1/24.
//

#ifndef EXPRESSION_TEMPLATES_H
#define EXPRESSION_TEMPLATES_H

#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <type_traits>
#include <stdexcept>


template<typename T>
struct Scalar {
private:
    const T& s;

public:
    constexpr explicit Scalar (const T& v) : s(v) {}

    constexpr const T& operator[] (const std::size_t) const {
        return s;
    }

    constexpr const T& operator() (const std::size_t, const std::size_t) const {
        return s;
    }

    [[nodiscard]] constexpr std::size_t size() const { return 0u; }
};

template<typename T>
struct ExprTraits {
    using ExprRef = const T&;
};

template<typename T>
struct ExprTraits<Scalar<T>> {
    using ExprRef = Scalar<T>;
};

/* ============================================
 * Expression Templates
 * =========================================== */


template <typename T, typename E1, typename E2>
class AddExp {
private:
    typename ExprTraits<E1>::ExprRef lhs;
    typename ExprTraits<E2>::ExprRef rhs;

public:
    AddExp(const E1& a, const E2& b) : lhs(a), rhs(b) {}

    T operator()(const std::size_t i, const std::size_t j) {
        return lhs(i, j) + rhs(i, j);
    }

    [[nodiscard]] std::size_t size() const {
        assert(lhs.size() == 0 or rhs.size() == 0 or lhs.size() == rhs.size());
        return lhs.size() != 0 ? lhs.size() : rhs.size();
    }
};




/* ============================================
 * Matrix Class
 * =========================================== */
template <typename T = double, std::size_t ROWS_ = 3, std::size_t COLS_ = 3>
class Matrix {
public:
    using value_t = T;

    Matrix() : m_data{} {}

    Matrix(Matrix const&) = default;

    ~Matrix() = default;

    Matrix& operator=(Matrix const& o) {
        if (&o != this) {
            copy(o);
        }
        return *this;
    }

    [[nodiscard]] static std::size_t size() {
        return ROWS * COLS;
    }

    T& operator()(const std::size_t i, const std::size_t j) { return m_data[j + (i * COLS_)]; }
    const T& operator()(const std::size_t i, const std::size_t j) const { return m_data[j + (i * COLS_)]; }

protected:
    void copy(const Matrix& o) {
        assert(ROWS == o.ROWS and COLS == o.COLS);
        for (std::size_t idx = 0; idx < size(); idx++) {
            m_data[idx] = o.m_data[idx];
        }

    }

public:
    static constexpr std::size_t ROWS = ROWS_;
    static constexpr std::size_t COLS = COLS_;
private:
    std::array<T, ROWS * COLS> m_data;
};


template<typename T>
Matrix<T> operator+(const Matrix<T>& a, const Matrix<T>& b) {
    assert(a.ROWS == b.ROWS and a.COLS == b.COLS);

}


#endif //EXPRESSION_TEMPLATES_H
