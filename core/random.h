//
// Created by cepheid on 7/31/24.
//

#ifndef RANDOM_H
#define RANDOM_H

#include <random>

// RNG class definition
template<template<typename> class Distribution, std::floating_point T>
class RandomNumberGenerator {
public:
  explicit RandomNumberGenerator(std::decay_t<T> low, std::decay_t<T> high) : distribution{low, high}, generator(init_mt_64()){}

  auto operator()() {
    return distribution(generator);
  }

private:
  Distribution<T> distribution;
  std::mt19937_64 generator;

  static inline std::mt19937_64 init_mt_64() {
    std::array<int, 624> seed_data{};
    std::random_device r; // @akis: this may be problematic for portability later, be careful
    std::generate_n(seed_data.data(), seed_data.size(), std::ref(r));
    std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
    std::mt19937_64 gen(seq);
    return gen;
  }
};

// using uniform_rng = RandomNumberGenerator<std::uniform_real_distribution, fptype, fptype>;
// using normal_rng = RandomNumberGenerator<std::normal_distribution, fptype, fptype>;


#endif //RANDOM_H
