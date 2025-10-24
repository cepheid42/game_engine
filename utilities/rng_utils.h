#ifndef TRIFORCE_RNG_UTILS_H
#define TRIFORCE_RNG_UTILS_H

#include <random>
#include <array>
#include <algorithm>

namespace tf::random {
  struct shuffleVector {
    std::random_device rd;
    std::mt19937_64 g;

    shuffleVector() : g(rd()) {}

    template<typename T>
    void operator() (std::vector<T> &vec, const size_t start, const size_t stop) {
      std::shuffle(vec.begin() + start, vec.begin() + stop, g);
    }
  };

  // RNG class definition
  template<typename T, template<typename> class Distribution, typename... Args>
  class RandomNumberGenerator {
  public:
    explicit RandomNumberGenerator(std::decay_t<Args>... args) : distribution{args...}, generator(init_mt_64()) {}
    auto operator()() { return distribution(generator); }

  private:
    Distribution<T> distribution;
    std::mt19937_64 generator;

    static std::mt19937_64 init_mt_64() {
      std::array<int, 624> seed_data{};
      std::random_device r; // @akis: this may be problematic for portability later, be careful
      std::generate_n(seed_data.data(), seed_data.size(), std::ref(r));
      std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
      const std::mt19937_64 gen(seq);
      return gen;
    }
    //
  }; // end class RandomNumberGenerator

  template <std::floating_point fp>
  using uniform_rng = RandomNumberGenerator<fp, std::uniform_real_distribution, fp, fp>;

  template<std::integral T>
  using uniform_int_rng = RandomNumberGenerator<T, std::uniform_int_distribution, T, T>;

  template <std::floating_point fp>
  using normal_rng = RandomNumberGenerator<fp, std::normal_distribution, fp, fp>;

} // end namespace tf::math::random

#endif //TRIFORCE_RNG_UTILS_H
