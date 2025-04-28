#ifndef SMALL_VECTOR_HPP
#define SMALL_VECTOR_HPP

#include <array>
#include <memory>
#include <cassert>

template<typename T, std::size_t N>
struct small_vector {
  using value_type = T;
  using difference_type = std::ptrdiff_t;
  using size_type = std::size_t;
  using reference = T&;
  using const_reference = const T&;
  using pointer = T*;
  using const_pointer = const T*;
  using iterator = T*;
  using const_iterator = const T*;

  std::array<T, N> data;
  std::size_t size_{};
  std::size_t cap = N;

  constexpr small_vector() = default;

  constexpr explicit small_vector(std::size_t n) : size_(n) {
    assert(n <= N);
    std::uninitialized_fill(std::begin(data), n, T{});
  }

  constexpr small_vector(const small_vector& other) : size_(other.size()) {
    static_assert(other.cap <= cap);
    size_ = other.size;
    copy_construct(std::begin(other), std::end(other), std::begin(data));
  }

  template<typename Iterator>
  constexpr small_vector(Iterator first, Iterator last) {
    size_ = std::distance(first, last);
    std::uninitialized_copy(first, last, get_data());
  }

  constexpr small_vector(small_vector&& other) noexcept : size_(other.size_) {
    *this = std::move(other);
  }

  constexpr small_vector(std::initializer_list<T> other) {
    assert(other.size() <= N);
    size_ = other.size();
    std::uninitialized_copy(std::begin(other), std::end(other), get_data());
  }

  constexpr ~small_vector() = default;

  constexpr small_vector& operator=(const small_vector& other) {
    if (this == &other) { return *this; }
    assign_copy(other);
    return *this;
  }

  constexpr small_vector& operator=(small_vector&& other) noexcept {
    if (this == &other) { return *this; }
    assign_move(other);
    return *this;
  }

  constexpr small_vector& operator=(std::initializer_list<T> other) {
    clear();
    size_ = other.size;
    copy_construct(std::begin(other), std::end(other), std::begin(data));
    return *this;
  }

  [[nodiscard]] constexpr T const* get_data() const { return data.data(); }
  constexpr T* get_data() { return data.data(); }

  [[nodiscard]] constexpr std::size_t size() const { return size_; }
  [[nodiscard]] constexpr bool empty() const { return size_ == 0;}

  constexpr iterator begin() { return get_data(); }
  [[nodiscard]] constexpr const_iterator begin() const { return get_data(); }

  constexpr iterator end() { return get_data() + size_; }
  [[nodiscard]] constexpr const_iterator end() const { return get_data() + size_; }

  constexpr iterator cbegin() { return get_data(); }
  [[nodiscard]] constexpr const_iterator cbegin() const { return get_data(); }

  constexpr iterator cend() { return get_data() + size_; }
  [[nodiscard]] constexpr const_iterator cend() const { return get_data() + size_; }

  constexpr void clear() {
    if (empty()) { return; }
    std::destroy(std::begin(data), std::end(data));
    size_ = 0;
  }

  template<typename Iterator>
  constexpr void construct_from_range(Iterator first, Iterator last) {
    std::uninitialized_copy(first, last, get_data());
  }

  template<typename... Args>
  constexpr T& emplace_back(Args&&... args) {
    auto* ptr = get_data() + size_;
    size_ += 1;
    std::construct_at(ptr, std::forward<Args>(args)...);
    return *ptr;
  }

  constexpr void push_back(const T& value) { emplace_back(value); }
  constexpr void push_back(T&& value) { emplace_back(std::move(value)); }

};

#endif //SMALL_VECTOR_HPP
