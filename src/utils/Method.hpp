#ifndef CDWPTTCHI2RECTOOL_UTILS_METHOD_HPP_
#define CDWPTTCHI2RECTOOL_UTILS_METHOD_HPP_

#include <type_traits>

#include "utils/RecPmtProp.hpp"

struct MethodTag {};

// These are "pure" tags, they are used for the template of Chi2, Predictor and CostFunction
struct FhtMethodTag : public MethodTag {};
struct TtMethodTag : public MethodTag {};

// These are "mixed" tags, they should only be used for the template of CostFunction
struct FhtTtMethodTag : public FhtMethodTag, public TtMethodTag {};

template<typename _Tag>
struct MethodTraits;

template<>
struct MethodTraits<FhtMethodTag> {
    typedef RecPmtProp value_type;
    typedef RecPmtTable vector_type;
    typedef RecPmtTable::iterator iterator;
    typedef RecPmtTable::const_iterator const_iterator;

    typedef double theo_value_type;
    typedef std::vector<double> theo_vector_type;
    typedef std::vector<double>::iterator theo_iterator;
    typedef std::vector<double>::const_iterator theo_const_iterator;
};

template<>
struct MethodTraits<TtMethodTag> {
    typedef vec3 value_type;
    typedef std::vector<vec3> vector_type;
    typedef std::vector<vec3>::iterator iterator;
    typedef std::vector<vec3>::const_iterator const_iterator;

    typedef vec3 theo_value_type;
    typedef std::vector<vec3> theo_vector_type;
    typedef std::vector<vec3>::iterator theo_iterator;
    typedef std::vector<vec3>::const_iterator theo_const_iterator;
};

#endif // CDWPTTCHI2RECTOOL_UTILS_METHOD_HPP_