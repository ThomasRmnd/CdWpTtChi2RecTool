#ifndef CDWPTTCHI2RECTOOL_PREDICTOR_PREDICTOR_HPP_
#define CDWPTTCHI2RECTOOL_PREDICTOR_PREDICTOR_HPP_

#include "utils/TrackParams.hpp"

#include "utils/Method.hpp"

#include <memory>

template<typename _Tag>
class Predictor : public IParamsHandler {

    static_assert(std::is_base_of<MethodTag, _Tag>::value, "Tag must derive from MethodTag");

public:

    typedef typename MethodTraits<_Tag>::const_iterator const_iterator;
    typedef typename MethodTraits<_Tag>::theo_iterator theo_iterator;

    virtual ~Predictor() = default;

    virtual void predict(const_iterator first, const_iterator last, theo_iterator theo, const double* params) = 0;

};

#endif // CDWPTTCHI2RECTOOL_PREDICTOR_PREDICTOR_HPP_