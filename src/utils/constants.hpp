#ifndef CDWPTTCHI2RECTOOL_UTILS_CONSTANTS_HPP_
#define CDWPTTCHI2RECTOOL_UTILS_CONSTANTS_HPP_

#include <cmath>

namespace constants {

// "pure" constants

template<typename _Tp>
inline constexpr _Tp pi_v = static_cast<_Tp>(M_PI);

inline constexpr double pi = pi_v<double>;

template<typename _Tp>
inline constexpr _Tp r_acrylic_v = static_cast<_Tp>(17700.0); // Radius of the acrylic (mm)

inline constexpr double r_acrylic = r_acrylic_v<double>;

template<typename _Tp>
inline constexpr _Tp r_cd_v = static_cast<_Tp>(20050.0); // Radius of the CD (mm)

inline constexpr double r_cd = r_cd_v<double>;

template<typename _Tp>
inline constexpr _Tp r_wp_v = static_cast<_Tp>(21750.0); // Radius of the WP (mm)

inline constexpr double r_wp = r_wp_v<double>;

template<typename _Tp>
inline constexpr _Tp h_wp_v = static_cast<_Tp>(44000.0); // Height of the WP (mm)

inline constexpr double h_wp = h_wp_v<double>;

template<typename _Tp>
inline constexpr _Tp c_v = static_cast<_Tp>(299.792458); // speed of light (mm/ns)

inline constexpr double c = c_v<double>;

template<typename _Tp>
inline constexpr _Tp n_ls_v = static_cast<_Tp>(1.485); // refractive index of the LS (1)

inline constexpr double n_ls = n_ls_v<double>;

template<typename _Tp>
inline constexpr _Tp n_w_v = static_cast<_Tp>(1.355); // refractive index of the water (1)

inline constexpr double n_w = n_w_v<double>;

// "mixed" constants

template<typename _Tp>
inline constexpr _Tp inv_c_v = static_cast<_Tp>(1.0) / c_v<_Tp>; // inverse of the speed of light (ns/mm)

inline constexpr double inv_c = inv_c_v<double>;

template<typename _Tp>
inline constexpr _Tp inv_c_ls_v = n_ls_v<_Tp> * inv_c_v<_Tp>; // inverse of the speed of light in the LS (ns/mm)

inline constexpr double inv_c_ls = inv_c_ls_v<double>;

template<typename _Tp>
inline constexpr _Tp inv_c_w_v = n_w_v<_Tp> * inv_c_v<_Tp>; // inverse of the speed of light in the water (ns/mm)

inline constexpr double inv_c_w = inv_c_w_v<double>;

template<typename _Tp>
inline _Tp tan_cherenkov_ls_v = std::sqrt(n_ls_v<_Tp> * n_ls_v<_Tp> - static_cast<_Tp>(1.0)); // tangent of the cherenkov angle in the LS (1)

inline double tan_cherenkov_ls = tan_cherenkov_ls_v<double>;

template<typename _Tp>
inline _Tp tan_cherenkov_w_v = std::sqrt(n_w_v<_Tp> * n_w_v<_Tp> - static_cast<_Tp>(1.0)); // tangent of the cherenkov angle in the water (1)

inline double tan_cherenkov_w = tan_cherenkov_w_v<double>;

template<typename _Tp>
inline _Tp inv_tan_cherenkov_ls_v = static_cast<_Tp>(1.0) / tan_cherenkov_ls_v<_Tp>; // inverse of the tangent of the cherenkov angle in the LS (1)

inline double inv_tan_cherenkov_ls = inv_tan_cherenkov_ls_v<double>;

template<typename _Tp>
inline _Tp inv_tan_cherenkov_w_v = static_cast<_Tp>(1.0) / tan_cherenkov_w_v<_Tp>; // inverse of the tangent of the cherenkov angle in the water (1)

inline double inv_tan_cherenkov_w = inv_tan_cherenkov_w_v<double>;

template<typename _Tp>
inline _Tp inv_c_tan_cherenkov_ls_v = inv_c_v<_Tp> * inv_tan_cherenkov_ls_v<_Tp>; // inverse of the tangent of the cherenkov angle in the LS multiplied by the speed of light (ns/mm)

inline double inv_c_tan_cherenkov_ls = inv_c_tan_cherenkov_ls_v<double>;

template<typename _Tp>
inline _Tp inv_c_tan_cherenkov_w_v = inv_c_v<_Tp> * inv_tan_cherenkov_w_v<_Tp>; // inverse of the tangent of the cherenkov angle in the water multiplied by the speed of light (ns/mm)

inline double inv_c_tan_cherenkov_w = inv_c_tan_cherenkov_w_v<double>;

template<typename _Tp>
inline constexpr _Tp n_ls_2_div_c_v = n_ls_v<_Tp> * n_ls_v<_Tp> * inv_c_v<_Tp>; // \frac{n_{ls}^2}{c} (ns/mm)

inline constexpr double n_ls_2_div_c = n_ls_2_div_c_v<double>;

template<typename _Tp>
inline constexpr _Tp n_w_2_div_c_v = n_w_v<_Tp> * n_w_v<_Tp> * inv_c_v<_Tp>; // \frac{n_w^2}{c} (ns/mm)

inline constexpr double n_w_2_div_c = n_w_2_div_c_v<double>;

template<typename _Tp>
inline constexpr _Tp one_minus_n_ls_2_div_c_v = (static_cast<_Tp>(1.0) - n_ls_v<_Tp> * n_ls_v<_Tp>) * inv_c_v<_Tp>; // \frac{1 - n_{ls}^2}{c} (ns/mm)

inline constexpr double one_minus_n_ls_2_div_c = one_minus_n_ls_2_div_c_v<double>;

template<typename _Tp>
inline constexpr _Tp one_minus_n_w_2_div_c_v = (static_cast<_Tp>(1.0) - n_w_v<_Tp> * n_w_v<_Tp>) * inv_c_v<_Tp>; // \frac{1 - n_w^2}{c} (ns/mm)

inline constexpr double one_minus_n_w_2_div_c = one_minus_n_w_2_div_c_v<double>;

} // namespace constants

#endif // CDWPTTCHI2RECTOOL_UTILS_CONSTANTS_HPP_