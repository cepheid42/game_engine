#ifndef TRAITS_HPP
#define TRAITS_HPP

#include <concepts>

// struct periodic_t {};
// struct pml_t {};
// struct null_t {};
//
// template<typename T> concept is_pml = std::derived_from<T, pml_t>;
// template<typename T> concept is_null = std::derived_from<T, null_t>;
// template<typename T> concept is_periodic = std::derived_from<T, periodic_t>;

enum class EMFace { X, Y, Z };
enum class EMSide { Lo, Hi };

template<EMFace F> concept is_XFace = F == EMFace::X;
template<EMFace F> concept is_YFace = F == EMFace::Y;
template<EMFace F> concept is_ZFace = F == EMFace::Z;
template<EMSide S> concept is_LoSide = S == EMSide::Lo;
template<EMSide S> concept is_HiSide = S == EMSide::Hi;

template<typename T> concept has_Ex = requires (T t) { t.Ex; t.Jx; };
template<typename T> concept has_Ey = requires (T t) { t.Ey; t.Jy; };
template<typename T> concept has_Ez = requires (T t) { t.Ez; t.Jz; };
template<typename T> concept has_Hx = requires (T t) { t.Hx; };
template<typename T> concept has_Hy = requires (T t) { t.Hy; };
template<typename T> concept has_Hz = requires (T t) { t.Hz; };

template<typename T> concept has_Ex_app = requires (T t) { t.Ex_app; t.Ex_total; };
template<typename T> concept has_Ey_app = requires (T t) { t.Ey_app; t.Ey_total; };
template<typename T> concept has_Ez_app = requires (T t) { t.Ez_app; t.Ez_total; };
template<typename T> concept has_Bx_app = requires (T t) { t.Bx_app; t.Hx_total; };
template<typename T> concept has_By_app = requires (T t) { t.By_app; t.Hy_total; };
template<typename T> concept has_Bz_app = requires (T t) { t.Bz_app; t.Hz_total; };

template<typename T> concept is_TMx = has_Ey<T> and has_Ez<T> and has_Hx<T>;
template<typename T> concept is_TMy = has_Ex<T> and has_Ez<T> and has_Hy<T>;
template<typename T> concept is_TMz = has_Ex<T> and has_Ey<T> and has_Hz<T>;

template<typename T> concept is_3D = is_TMx<T> and is_TMy<T> and is_TMz<T>;
template<typename T> concept is_2D = (is_TMx<T> or is_TMy<T> or is_TMz<T>) and !is_3D<T>;
// template<typename T> concept is_1D = (has_Ex<T> and (has_Hy<T> or has_Hz<T>)) or
//                                      (has_Ey<T> and (has_Hx<T> or has_Hz<T>)) or
//                                      (has_Ez<T> and (has_Hx<T> or has_Hy<T>));

#endif //TRAITS_HPP
