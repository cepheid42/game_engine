//
// Created by cepheid on 7/11/24.
//

#ifndef ITERATOR_H
#define ITERATOR_H

#include <cstddef>  // for std::ptrdiff_t
#include <iterator> // for std::foward_iterator_tag

template<typename T>
class MultiIterator {
public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = T;
    using difference_type = std::ptrdiff_t;
    using pointer = T*;
    using reference = T&;

    MultiIterator() = default;
    MultiIterator(pointer f,
                  pointer d1,
                  pointer d2,
                  pointer j,
                  difference_type f_str,
                  difference_type d1_str,
                  difference_type d2_str,
                  difference_type j_str)
    : f_ptr(f),
      d1_ptr(d1),
      d2_ptr(d2),
      j_ptr(j),
      f_stride(f_str),
      d1_stride(d1_str),
      d2_stride(d2_str),
      j_stride(j_str)
    {}

    struct Result
    {

    };

    MultiIterator& operator++() {
        f_ptr += f_stride;
        d1_ptr += d1_stride;
        d2_ptr += d2_stride;
        j_ptr += j_stride;
        return *this;
    }

    MultiIterator operator++(int) {
        MultiIterator temp = *this;
        f_ptr += f_stride;
        d1_ptr += d1_stride;
        d2_ptr += d2_stride;
        j_ptr += j_stride;
        return temp;
    }

    friend bool operator==(const MultiIterator& lhs, const MultiIterator& rhs) {
        const auto same_f = lhs.f_ptr == rhs.f_ptr;
        const auto same_d1 = lhs.d1_ptr == rhs.d1_ptr;
        const auto same_d2 = lhs.d2_ptr == rhs.d2_ptr;
        const auto same_j = lhs.j_ptr == rhs.j_ptr;
        return same_f and same_d1 and same_d2 and same_j;
    }

    friend bool operator!=(const MultiIterator& lhs, const MultiIterator& rhs) {
        return !(lhs == rhs);
    }

private:
    pointer f_ptr;
    pointer d1_ptr;
    pointer d2_ptr;
    pointer j_ptr;

    difference_type f_stride;
    difference_type d1_stride;
    difference_type d2_stride;
    difference_type j_stride;
};

#endif //ITERATOR_H
