#ifndef CDWPTTCHI2RECTOOL_UTILS_CHI2_HPP_
#define CDWPTTCHI2RECTOOL_UTILS_CHI2_HPP_

#include <numeric>

#include "utils/RecPmtProp.hpp"

inline double chi2FhtMethod(RecPmtTable::const_iterator first, RecPmtTable::const_iterator last, std::vector<double>::const_iterator theo) {
    return std::transform_reduce(first, last, theo, 0.0, std::plus<double>(), [](const RecPmtProp& pmt, double value) {
        return std::pow((pmt.fht - value) * pmt.inv_res, 2.0);
    });
}

inline double chi2TtMethod(std::vector<vec3>::const_iterator first, std::vector<vec3>::const_iterator last, std::vector<vec3>::const_iterator theo) {
    return std::transform_reduce(first, last, theo, 0.0, std::plus<double>(), [](const vec3& hit, const vec3& value) {
        return mag2(hit - value);
    });
}

#endif // CDWPTTCHI2RECTOOL_UTILS_CHI2_HPP_