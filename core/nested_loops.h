//
// Created by cepheid on 7/12/24.
//

#ifndef NESTED_LOOPS_H
#define NESTED_LOOPS_H

#include <array>

// template<size_t DEPTH>
// struct NestedLoopGen {
//     using iterator_category = std::forward_iterator_tag;
//     using value_type = size_t;
//     using difference_type = std::ptrdiff_t;
//     using pointer = size_t*;
//     using reference = size_t&;
//
//     NestedLoopGen() = default;
//
//     NestedLoopGen(const std::array<size_t, DEPTH>& min, const std::array<size_t, DEPTH>& max)
//     : index(DEPTH - 1),
//       min(min),
//       max(max),
//       slots(min)
//     {}
//
//     reference operator*() const { return slots; }
//
//     NestedLoopGen& operator++() {
//         // Increment
//         ++slots[DEPTH - 1];
//
//         // Carry
//         while (slots[index] == max[index]) {
//             // Overflow, we're done
//             if (index == 0) { return; }
//
//             slots[index] = min[index];
//             ++slots[--index];
//         }
//         index = DEPTH - 1;
//     }
//
//     size_t index;
//     std::array<size_t, DEPTH> min;
//     std::array<size_t, DEPTH> max;
//     std::array<size_t, DEPTH> slots;
// };

template<size_t DEPTH, typename Func>
void IterativeNestedLoop(const std::array<size_t, DEPTH>& min, const std::array<size_t, DEPTH>& max, Func action) {
    // Initialize the slots to hold the current iteration value for each depth
    auto slots = min;

    size_t index = DEPTH - 1;
    // while (true) {
        action(slots, DEPTH);

        // Increment
        ++slots[DEPTH - 1];

        // Carry
        while (slots[index] == max[index]) {
            // Overflow, we're done
            if (index == 0) { return; }

            slots[index] = min[index];
            ++slots[--index];
        }
        index = DEPTH - 1;
    // }
}

#endif //NESTED_LOOPS_H
