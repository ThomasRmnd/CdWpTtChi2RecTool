#ifndef CDWPTTCHI2RECTOOL_PREDICTOR_PREDICTOR_HPP_
#define CDWPTTCHI2RECTOOL_PREDICTOR_PREDICTOR_HPP_

#include "utils/Configurable.hpp"
#include "utils/Method.hpp"
#include "utils/TrackParams.hpp"

/**
 * @class Predictor
 * @tparam _MethodTag Method tag
 * 
 * @brief Base class to calculate expected data
 */
template<typename _MethodTag>
class Predictor : public IParamsHandler {

    static_assert(std::is_base_of<MethodTag, _MethodTag>::value, "Tag must derive from MethodTag");

public:

    typedef typename MethodTraits<_MethodTag>::const_iterator const_iterator;
    typedef typename MethodTraits<_MethodTag>::theo_iterator theo_iterator;

    virtual ~Predictor() = default;

    /**
     * @brief Calculate the expected data
     * 
     * @param first the beginning iterator of the experimental data 
     * @param last the end iterator of the experimental data
     * @param theo the beginning iterator of the expected data
     * @param params track parameters
     */
    virtual void predict(const_iterator first, const_iterator last, theo_iterator theo, const double* params) = 0;

};

#endif // CDWPTTCHI2RECTOOL_PREDICTOR_PREDICTOR_HPP_