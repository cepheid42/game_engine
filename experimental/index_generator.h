//
// Created by cepheid on 7/12/24.
//

#ifndef INDEX_GENERATOR_H
#define INDEX_GENERATOR_H

#include <fstream>
#include <cmath>
#include <iostream>
#include <sstream>
#include <functional>

// template<size_t DEPTH, typename Func>
// void IterativeNestedLoop(const std::array<size_t, DEPTH>& min, const std::array<size_t, DEPTH>& max, Func action) {
//     // Initialize the slots to hold the current iteration value for each depth
//     auto slots = min;
//
//     while (true) {
//         action(slots, DEPTH);
//
//         // Increment
//         ++slots[DEPTH - 1];
//
//         // Carry
//         size_t index = DEPTH - 1;
//         while (slots[index] == max[index]) {
//             // Overflow, we're done
//             if (index == 0) { return; }
//
//             slots[index] = min[index];
//             ++slots[--index];
//         }
//     }
// }

template<size_t DEPTH>
struct NestedLoop
{
    NestedLoop(const std::array<size_t, DEPTH>& min, const std::array<size_t, DEPTH>& max)
    : index(DEPTH - 1),
      count(0),
      max_iter(1),
      min(min),
      max(max),
      slots(min)
    {
        for (size_t i = 0; i < DEPTH; i++) {
            const auto val = max[i] - min[i];
            max_iter *= val;
        }
    }

    auto operator()() {
        auto result = slots;

        // Increment
        ++slots[DEPTH - 1];

        // Carry
        index = DEPTH - 1;
        while (slots[index] == max[index]) {
            // Overflow, we're done
            if (index == 0) { break; }

            slots[index] = min[index];
            ++slots[--index];
        }

        count++;

        return result;
    }

    void reset()
    {
        count = 0;
        slots = min;
    }

    explicit operator bool() {
        if (count != max_iter) { return true; }

        reset();
        return false;
    }

    size_t index;
    size_t count;
    size_t max_iter;
    std::array<size_t, DEPTH> min;
    std::array<size_t, DEPTH> max;
    std::array<size_t, DEPTH> slots;
};



void DoStuff(const std::array<size_t, 3>& slots, const size_t depth) {
    std::stringstream values;
    for(size_t i = 0; i < depth; i++)
    {
        values << slots[i] << ", ";
    }
    std::cout << values.str() << std::endl;
}

void repeat(NestedLoop<3>& looper) {
    while (looper) {
        auto results = looper();
        DoStuff(results, 3);
    }
}

/*
int main() {
    std::array<size_t, 3> min = {1, 1, 1};
    std::array<size_t, 3> max = {5, 5, 5};

    NestedLoop looper(min, max);

    for (size_t i = 0; i < 3; i++) {
        std::cout << "New Looper\n";
        repeat(looper);
        std::cout << std::endl;
    }
}
*/
#endif //INDEX_GENERATOR_H
