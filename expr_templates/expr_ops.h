//
// Created by cepheid on 8/19/24.
//

#ifndef EXPR_OPS_H
#define EXPR_OPS_H

#include <cassert>
#include <cstddef>

#include "expr_array.h"

template<typename T>
class A_Scalar;

// primary template
template<typename T>
struct A_Traits {
  using ExprRef = const T&;
};

// partial specialization for scalars
template<typename T>
struct A_Traits<A_Scalar<T>> {
  using ExprRef = A_Scalar<T>;
};

// template<typename T> struct expr_holder { using type = T; };
// template<typename T> struct expr_holder<const T> : expr_holder<T> {};
// template<typename T> struct expr_holder<T&> { using type = const T&; };
// template<typename T> struct expr_holder<T&&> { using type = T; };
//
// template<typename T> using expr_holder_t = typename expr_holder<T>::type;

// class for objects that represent the addition of two operands
template<typename T, typename OP1, typename OP2>
class A_Add {
  // const OP1& op1; // first operand
  // const OP2& op2; // second operand

  typename A_Traits<OP1>::ExprRef op1; // first operand
  typename A_Traits<OP2>::ExprRef op2; // second operand

public:
  // constructor initializes references to operands
  A_Add(const OP1& a, const OP2& b) : op1(a), op2(b) {
  }

  // compute sum when value requested
  T operator[](const std::size_t idx) const { return op1[idx] + op2[idx]; }

  // size is maximum size
  [[nodiscard]] std::size_t size() const {
    assert(op1.size()==0 || op2.size()==0 || op1.size()==op2.size());
    return op1.size() != 0 ? op1.size() : op2.size();
  }
};

// class for objects that represent the addition of two operands
template<typename T, typename OP1, typename OP2>
class A_Sub {
  // const OP1& op1; // first operand
  // const OP2& op2; // second operand

  typename A_Traits<OP1>::ExprRef op1; // first operand
  typename A_Traits<OP2>::ExprRef op2; // second operand

public:
  // constructor initializes references to operands
  A_Sub(const OP1& a, const OP2& b) : op1(a), op2(b) {
  }

  // compute sum when value requested
  T operator[](const std::size_t idx) const { return op1[idx] - op2[idx]; }

  // size is maximum size
  [[nodiscard]] std::size_t size() const {
    assert(op1.size()==0 || op2.size()==0 || op1.size()==op2.size());
    return op1.size() != 0 ? op1.size() : op2.size();
  }
};

// class for objects that represent the multiplication of two operands
template<typename T, typename OP1, typename OP2>
class A_Mult {
  // const OP1& op1; // first operand
  // const OP2& op2; // second operand

  typename A_Traits<OP1>::ExprRef op1; // first operand
  typename A_Traits<OP2>::ExprRef op2; // second operand

public:
  // constructor initializes references to operands
  A_Mult(const OP1& a, const OP2& b) : op1(a), op2(b) {
  }

  // compute product when value requested
  T operator[](const std::size_t idx) const { return op1[idx] * op2[idx]; }

  // size is maximum size
  [[nodiscard]] std::size_t size() const {
    assert(op1.size()==0 || op2.size()==0 || op1.size()==op2.size());
    return op1.size() != 0 ? op1.size() : op2.size();
  }
};

template<typename T, typename A1, typename A2>
class A_Subscript {
  A1& a1; // first operand
  const A2& a2; // second operand

public:
  A_Subscript(A1& a, const A2& b) : a1(a), a2(b) {}

  T& operator[](const std::size_t idx) { return a1[a2[idx]]; }
  const T& operator[](const std::size_t idx) const { return a1[a2[idx]]; }

  [[nodiscard]] std::size_t size() const { return a2.size(); }
};

template<typename T, typename R>
template<typename T2, typename R2>
decltype(auto)
Array<T, R>::operator[](const Array<T2, R2>& b) const {
  return Array<T, A_Subscript<T, R, R2>>(A_Subscript<T, R, R2>((*this).rep(), b.rep()));
}

template<typename T, typename R>
template<typename T2, typename R2>
decltype(auto)
Array<T, R>::operator[](const Array<T2, R2>& b) {
  return Array<T, A_Subscript<T, R, R2>>(A_Subscript<T, R, R2>((*this).rep(), b.rep()));
}

// addition of two Arrays:
template<typename T, typename R1, typename R2>
Array<T, A_Add<T, R1, R2>>
operator+(const Array<T, R1>& a, const Array<T, R2>& b) {
  return Array<T, A_Add<T, R1, R2>>(A_Add<T, R1, R2>(a.rep(), b.rep()));
}

// addition of scalar and Array:
template<typename T, typename R2>
Array<T, A_Add<T, A_Scalar<T>, R2>>
operator+(const T& s, const Array<T, R2>& b) {
  return Array<T, A_Add<T, A_Scalar<T>, R2>>
    (A_Add<T, A_Scalar<T>, R2>(A_Scalar<T>(s), b.rep()));
}

template<typename T, typename R2>
Array<T, A_Add<T, A_Scalar<T>, R2>>
operator+(const Array<T, R2>& b, const T& s) {
  return s + b;
}

// subtraction of two Arrays:
template<typename T, typename R1, typename R2>
Array<T, A_Sub<T, R1, R2>>
operator-(const Array<T, R1>& a, const Array<T, R2>& b) {
  return Array<T, A_Sub<T, R1, R2>>(A_Sub<T, R1, R2>(a.rep(), b.rep()));
}

// subtraction of array and scalar
template<typename T, typename R2>
Array<T, A_Sub<T, A_Scalar<T>, R2>>
operator-(const T& s, const Array<T, R2>& b) {
  return Array<T, A_Sub<T, A_Scalar<T>, R2>>(A_Sub<T, A_Scalar<T>, R2>(A_Scalar<T>(s), b.rep()));
}

template<typename T, typename R2>
Array<T, A_Sub<T, R2, A_Scalar<T>>>
operator-(const Array<T, R2>& b, const T& s) {
  return Array<T, A_Sub<T, R2, A_Scalar<T>>>(A_Sub<T, R2, A_Scalar<T>>(b.rep(), A_Scalar<T>(s)));
}

// multiplication of two Arrays:
template<typename T, typename R1, typename R2>
Array<T, A_Mult<T, R1, R2>>
operator*(const Array<T, R1>& a, const Array<T, R2>& b) {
  return Array<T, A_Mult<T, R1, R2>>(A_Mult<T, R1, R2>(a.rep(), b.rep()));
}

// multiplication of scalar and Array:
template<typename T, typename R2>
Array<T, A_Mult<T, A_Scalar<T>, R2>>
operator*(const T& s, const Array<T, R2>& b) {
  return Array<T, A_Mult<T, A_Scalar<T>, R2>>
    (A_Mult<T, A_Scalar<T>, R2>(A_Scalar<T>(s), b.rep()));
}

template<typename T, typename R2>
Array<T, A_Mult<T, A_Scalar<T>, R2>>
operator*(const Array<T, R2>& b, const T& s) {
  return s * b;
}

#endif //EXPR_OPS_H
