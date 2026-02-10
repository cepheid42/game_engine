#ifndef C_STRINGS_H
#define C_STRINGS_H

#include <array>
#include <string>

using namespace std::literals;

// https://fekir.info/post/compile-time-string-operations/
// Concatenation Operations
template<typename T, std::size_t N1, std::size_t N2, std::size_t... P1, std::size_t... P2>
constexpr auto concat_impl(const std::array<T, N1>& arr1, const std::array<T, N2>& arr2,
                           std::index_sequence<P1...>, std::index_sequence<P2...>) {
   return std::array<T, N1 + N2>({arr1[P1]..., arr2[P2]...});
}

template<typename T, std::size_t N1, std::size_t N2>
constexpr auto concat(const std::array<T, N1>& arr1, const std::array<T, N2>& arr2) {
   return concat_impl(arr1, arr2, std::make_index_sequence<N1>{}, std::make_index_sequence<N2>{});
}

static_assert(concat(std::array{1,2}, std::array{3,4}).size()==4,"");
static_assert(concat(concat(std::array{1,2}, std::array{1}), std::array{3,4}).size()==5,"");


// Trim Operations
template<std::size_t O, typename T, std::size_t N, std::size_t... P>
constexpr auto trim_right_impl(const std::array<T, N>& arr, std::index_sequence<P...>) {
   return std::array<T, N - O>({arr[P]...});
}

template<std::size_t O, typename T, std::size_t N>
constexpr auto trim_right(const std::array<T, N>& arr) {
   return trim_right_impl<O>(arr, std::make_index_sequence<N - O>{});
}

static_assert(trim_right<1>(std::to_array({1,2,3})).size() == 2);
static_assert(trim_right<1>(std::to_array({1,2,3}))[0] == 1);
static_assert(trim_right<1>(std::to_array({1,2,3}))[1] == 2);

// c_string class
template<std::size_t N>
struct c_string {
   using array_t = std::array<char, N + 1>;

   template<std::size_t N2>
   friend struct c_string;

   std::array<char, N + 1> data;

   template<std::size_t... I>
   constexpr c_string(const array_t& arr, std::index_sequence<I...>) noexcept : data{ arr[I]..., '\0' } {}
   constexpr c_string(const array_t& arr) noexcept : c_string(arr, std::make_index_sequence<N>{}) {}
   constexpr c_string(const char (&arr)[N + 1]) noexcept : c_string(std::to_array(arr), std::make_index_sequence<N>{}) {}

   template <std::size_t N2>
   constexpr auto operator+(const c_string<N2>& s2) const noexcept -> c_string<N + N2>{
      return c_string<N + N2>(concat(trim_right<1>(data), s2.data));
   }

   template<std::size_t N2>
   constexpr friend auto operator+(const c_string<N>& s1, const char (&s2)[N2]) noexcept -> c_string<N + N2 - 1> {
      return s1 + c_string<N2 - 1>{s2};
   }

   template<std::size_t N2>
   constexpr friend auto operator+(const char (&s1)[N2], const c_string<N>& s2) noexcept -> c_string<N + N2 - 1> {
      return c_string<N2 - 1>{s1} + s2;
   }

   static constexpr std::size_t size() noexcept { return N; }
   constexpr char operator[](int i) const { return data[i]; }

   constexpr operator std::string_view() { return std::string_view{data.data()}; }
   constexpr const char* c_str() const noexcept { return data.data(); }
};

// factory function
template<int NP1>
constexpr auto literal(const char (&lit)[NP1]) noexcept -> c_string<NP1 - 1> { return c_string<NP1 - 1>(lit); }

static_assert((literal("1")+literal("23")).size() == 3, "");
static_assert((literal("")+literal("23")).size() == 2, "");
static_assert((literal("")+literal("")).size() == 0, "");
static_assert((literal("1")+literal("23"))[0] == '1', "");
static_assert((literal("1")+literal("23"))[1] == '2', "");
static_assert((literal("1")+literal("23"))[2] == '3', "");
static_assert((literal("1")+"23").size() == 3, "");
static_assert((literal("1")+"").size() == 1, "");
static_assert(("1" + literal("23")).size() == 3, "");
static_assert(("" + literal("23")).size() == 2, "");

// Conversion functions
constexpr int num_digits(int x) {
   return x < 0  ? 1 + num_digits(-x) :
          x < 10 ? 1 :
                   1 + num_digits(x / 10);
}

static_assert(num_digits(-1) == 2, "");
static_assert(num_digits(-0) == 1, "");
static_assert(num_digits( 0) == 1, "");
static_assert(num_digits(10) == 2, "");

constexpr char digit_to_char(int c) { return c + '0'; }

template<int value>
constexpr auto to_string_1(std::true_type) -> c_string<1> {
   return c_string<1>({digit_to_char(value), '\0'});
}

template<int value>
constexpr auto to_string_1(std::false_type) -> c_string<2> {
   return c_string<2>({'-', digit_to_char(-value), '\0'});
}

template<int value>
constexpr auto to_string_i(std::true_type) -> c_string<num_digits(value)> {
   return to_string_1<value>(std::integral_constant<bool, (value >= 0)>());
}

template<int value> constexpr auto to_string() -> c_string<num_digits(value)>;

template<int value>
constexpr auto to_string_i(std::false_type) -> c_string<num_digits(value)> {
   return to_string<value / 10>() + to_string_1<std::abs(value) % 10>(std::true_type{});
}

template<int value>
constexpr auto to_string() -> c_string<num_digits(value)> {
   return to_string_i<value>(std::integral_constant<bool, num_digits(value) == 1 or (num_digits(value) == 2 and not (value > 0))>());
}

// static_assert(to_string<-9999999>() == "-9999999"sv, "");
// static_assert(to_string<123>() == "123"sv, "");


template<bool> constexpr auto getTypeStr() { return literal("b"); }
template<typename> constexpr auto getTypeStr() { return literal("?"); }
template<std::floating_point> constexpr auto getTypeStr() { return literal("f"); }
template<std::signed_integral> constexpr auto getTypeStr() { return literal("i"); }
template<std::unsigned_integral> constexpr auto getTypeStr() { return literal("u"); }

#endif //C_STRINGS_H